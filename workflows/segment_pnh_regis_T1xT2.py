#!/usr/bin/env python3
"""
    PNH anatomical segmentation pipeline given by Regis Trapeau wrapped in
    Nipype.
    
    Description
    --------------
    TODO :/
    
    Arguments
    -----------
    -data:
        Path to the BIDS directory that contain subjects' MRI data.
        
    -out:
        Nipype's processing directory. 
        It's where all the outputs will be saved.
        
    -subjects:
        IDs list of subjects to process.

    -ses
        session (leave blank if None)
    
    Example
    ---------
    python segment_pnh_regis_T1xT2.py -data /home/bastien/work/data/local_primavoice -out /home/bastien/work/projects/macapype/local_tests -template /home/bastien/work/data/templates/inia19/inia19-t1-brain.nii -priors /home/bastien/work/data/templates/inia19/inia19-prob_0.nii /home/bastien/work/data/templates/inia19/inia19-prob_1.nii /home/bastien/work/data/templates/inia19/inia19-prob_2.nii -ses 01 -sub Maga
    
    Resources files
    -----------------
    Brain templates are required to run this segmentation pipeline. Please,
    download the resources and specify the unzipped directory with the
    '-resources' argument in the command line.
    Resources are here:
    # TODO: find a permanent place
    https://cloud.int.univ-amu.fr/index.php/s/8bCJ5CWWPfHRyHs
    
    Requirements
    --------------
    This workflow use:
        - FSL
    
"""

# Authors : David Meunier (david.meunier@univ-amu.fr)
#           Bastien Cagna (bastien.cagna@univ-amu.fr)
#           Regis Trapeau (regis.trapeau@univ-amu.fr)

import nilearn as ni
import nibabel as nb
import os
import os.path as op
import argparse

from bids.layout import BIDSLayout

import nipype.interfaces.io as nio
import nipype.interfaces.utility as niu
import nipype.pipeline.engine as pe
import nipype.interfaces.ants as ants
import nipype.interfaces.fsl as fsl

from macapype.nodes.preproc import average_align
from macapype.nodes.bash_regis import T1xT2BET, T1xT2BiasFieldCorrection, \
                                      IterREGBET
from macapype.pipelines.extract_brain import create_old_segment_extraction_pipe
from macapype.utils.utils_tests import load_test_data


fsl.FSLCommand.set_default_output_type('NIFTI_GZ')


def gunzip(filename):
    import subprocess

    if filename[-3:] == ".gz":
        subprocess.check_output("gunzip " + filename, shell=True)
    else:
        ValueError("Non GZip file given")
    return filename[:-3]


def format_spm_priors(priors, fname="merged_tissue_priors.nii",
                      directory=None):
    """
    Arguments
    =========
    priors: str or list

    fname: str
        Filename of the concatenated 4D Nifti image

    directory: str or None
        If None, the directory of the first file listed in prios is used.
    """
    if isinstance(priors, str):
        img = nb.load(priors)
        if len(img.shape) == 4 and 3 <= img.shape[3] <= 6:
            return priors
        else:
            raise ValueError(
                "Given Nifti is 3D while 4D expected or do not have between 3 "
                "and 6 maps."
            )
    elif isinstance(priors, list):
        imgs = []
        for f in priors:
            if directory is None:
                directory = op.split(f)[0]
            imgs.append(nb.load(f))
        fmt_image = ni.image.concat_imgs(imgs)

        new_img_f = op.join(directory, fname)
        print(new_img_f)
        nb.save(fmt_image, new_img_f)
        return new_img_f
    raise ValueError(
        "Priors must be one or a list of paths to a Nifti images"
    )


def create_infosource(subject_ids):
    infosource = pe.Node(
        interface=niu.IdentityInterface(fields=['subject_id']),
        name="infosource"
    )
    infosource.iterables = [('subject_id', subject_ids)]

    return infosource


def create_datasource(data_dir, subjects=[], sessions=[], acqs=[]):
    """ Create a datasource node that have iterables following BIDS format """
    bids_datasource = pe.Node(
        interface=nio.BIDSDataGrabber(),
        name='bids_datasource'
    )

    bids_datasource.inputs.base_dir = data_dir
    bids_datasource.inputs.output_query = {
        'T1': {
            "datatype": "anat", "suffix": "T1w",
            "extensions": ["nii", ".nii.gz"]
        },
        'T2': {
            "datatype": "anat", "suffix": "T2w",
            "extensions": ["nii", ".nii.gz"]
        }
    }

    layout = BIDSLayout(data_dir)

    # Verbose
    print("BIDS layaout:", layout)
    print("\t", layout.get_subjects())
    print("\t", layout.get_sessions())
    print("\t", layout.get_acquisitions())

    subjects = layout.get_subjects() if len(subjects) == 0 else subjects
    acqs = layout.get_acquistions() if acqs and len(acqs) == 0 else acqs
    sessions = layout.get_sessions() if len(sessions) == 0 else sessions

    iterables = []
    iterables.append(('subject', subjects))
    iterables.append(('session', sessions))
    if acqs:
        iterables.append(('acquisition', acqs))
    #    # Add iteration over subjects
    #    if len(layout.get_subjects()) == 1:
    #         bids_datasource.inputs.subject = layout.get_subjects()[0]
    #    else:
    #        iterables.append(('subject', layout.get_subjects()))

    #    #  Add iteration over sessions
    #    if len(layout.get_sessions()) == 1:
    #         bids_datasource.inputs.session = layout.get_sessions()[0]
    #    else:
    #        iterables.append(('session', layout.get_sessions()))

    #    #  Add iteration over acquisitions
    #    if len(layout.get_acquisitions()) == 1:
    #         bids_datasource.inputs.acquisitions = layout.get_acquisitions()[0]
    #    else:
    #        iterables.append(('acquisition', layout.get_acquisitions()))

    #    if len(iterables) > 0:
    bids_datasource.iterables = iterables

    return bids_datasource
# def create_datasource(data_dir, sess):
#     datasource = pe.Node(
#         interface=nio.DataGrabber(infields=['subject_id'], outfields=['T1','T2']),
#         name='datasource'
#     )
#     datasource.inputs.base_directory = data_dir
#     datasource.inputs.template = 'sub-%s/{}/anat/sub-%s_{}_%s.nii'.format(sess, sess)
#     datasource.inputs.template_args = dict(
#         T1=[['subject_id', 'subject_id', "T1w"]],
#         T2=[['subject_id', 'subject_id', "T2w"]],
#     )
#     datasource.inputs.sort_filelist = True
#
#     return datasource


def create_bids_datasource(data_dir):

    bids_datasource = pe.Node(
        interface=nio.BIDSDataGrabber(),
        name='bids_datasource'
    )

    bids_datasource.inputs.base_dir = data_dir

    bids_datasource.inputs.output_query = {
        'T1': {
            "datatype": "anat",
            "suffix": "T1w",
            "extensions": ["nii", ".nii.gz"]
        },
        'T2': {
            "datatype": "anat",
            "suffix": "T2w",
            "extensions": ["nii", ".nii.gz"]
        }
    }

    layout = BIDSLayout(data_dir)
    print(layout)
    print(layout.get_subjects())
    print(layout.get_sessions())

    iterables = []
    if len(layout.get_subjects()) == 1:
         bids_datasource.inputs.subject = layout.get_subjects()[0]
    else:
        iterables.append(('subject', layout.get_subjects()[:2]))

    if len(layout.get_sessions()) == 1:
         bids_datasource.inputs.session = layout.get_sessions()[0]
    else:
        iterables.append(('session', layout.get_sessions()[:2]))

    if len(iterables):
        bids_datasource.iterables = iterables

    return bids_datasource


def create_segment_pnh_T1xT2(brain_template, priors,
                             name='T1xT2_segmentation_pipeline'):
    print(brain_template)
    print(priors)
    print("node name: ", name)

    # Creating pipeline
    seg_pipe = pe.Workflow(name=name)

    # Creating input node
    inputnode = pe.Node(
        niu.IdentityInterface(fields=['T1','T2', 'priors']),
        name='inputnode'
    )

    # Brain extraction (unused) + Cropping
    bet = pe.Node(T1xT2BET(m=True, aT2=True, c=10), name='bet')
    bet.n = 2
    seg_pipe.connect(inputnode, 'T1', bet, 't1_file')
    seg_pipe.connect(inputnode, 'T2', bet, 't2_file')

    # Bias correction of cropped images
    debias = pe.Node(T1xT2BiasFieldCorrection(), name='debias')
    seg_pipe.connect(bet, 't1_cropped_file', debias, 't1_file')
    seg_pipe.connect(bet, 't2_cropped_file', debias, 't2_file')
    seg_pipe.connect(bet, 'mask_file', debias, 'b')

    # Iterative registration to the INIA19 template
    reg = pe.Node(IterREGBET(), name='reg')
    reg.inputs.refb_file = brain_template
    seg_pipe.connect(debias, 't1_debiased_file', reg, 'inw_file')
    seg_pipe.connect(debias, 't1_debiased_brain_file', reg, 'inb_file')

    # Compute brain mask using old_segment of SPM and postprocessing on
    # tissues' masks
    extract_brain = create_old_segment_extraction_pipe(priors)
    seg_pipe.connect(reg, ('warp_file', gunzip), extract_brain, 'inputnode.T1')

    return seg_pipe


###############################################################################
def create_main_workflow(data_dir, process_dir, subject_ids, sessions,
                         acquisitions, template, priors):
    """ """
    main_workflow = pe.Workflow(name="T1xT2_processing_workflow")
    main_workflow.base_dir = process_dir

    # # Infosource
    # if subject_ids is None or sessions is None:
    #     print('adding BIDS data source')
    #     datasource = create_bids_datasource(data_dir)
    # else:
    #     print('adding info source and data source')
    #     infosource = create_infosource(subject_ids)
    #     # Data source
    #     datasource = create_datasource(data_dir, sessions)
    #     # connect
    #     main_workflow.connect(
    #         infosource, 'subject_id', datasource, 'subject_id')
    datasource = create_datasource(data_dir, subject_ids, sessions,
                                   acquisitions)

    segment_pnh = create_segment_pnh_T1xT2(template, priors)

    main_workflow.connect(
        datasource, ('T1', average_align), segment_pnh, 'inputnode.T1')
    main_workflow.connect(
        datasource, ('T2', average_align), segment_pnh, 'inputnode.T2')

    return main_workflow


################################################################################
def main(data_path, main_path, subjects, sessions, acquisitions,
         template, priors):
    data_path = op.abspath(data_path)

    if not op.isdir(main_path):
        os.makedirs(main_path)

    # main_workflow
    print("Initialising the pipeline...")
    wf = create_main_workflow(
        data_dir=data_path,
        process_dir=main_path,
        subject_ids=subjects,
        sessions=sessions,
        acquisitions=acquisitions,
        template=template,
        priors=format_spm_priors(priors, directory=main_path)
    )
    # wf.write_graph(graph2use="colored")
    wf.config['execution'] = {'remove_unnecessary_outputs': 'false'}
    print('The PNH segmentation pipeline is ready')
    
    print("Start to process")
    wf.run()
    # wf.run(plugin='MultiProc', plugin_args={'n_procs' : 2})


if __name__ == '__main__':
    # Command line parser
    parser = argparse.ArgumentParser(
        description="PNH segmentation pipeline from Regis Trapeau")
    parser.add_argument("-data", dest="data", type=str, required=True,
                        help="Directory containing MRI data (BIDS)")
    parser.add_argument("-out", dest="out", type=str,
                        help="Output directory", required=True)
    parser.add_argument("-ses", dest="ses", type=str, nargs='+',
                        help="Sessions ID")
    parser.add_argument("-acq", dest="acq", type=str, nargs='+',
                        help="Acquisitions ID")
    parser.add_argument("-subjects", dest="subjects", type=str, nargs='+',
                        help="Subjects' ID", required=False)
    parser.add_argument("-template", dest="template", type=str,
                        default=None, help="Anatomical template")
    parser.add_argument("-priors", dest="priors", type=str, nargs='+',
                        default=None, help="Tissues probability maps")
    args = parser.parse_args()

    if args.template is None and args.priors is None:
        inia_dir = load_test_data("inia19")
        args.template = op.join(inia_dir, "inia19-t1-brain.nii")
        args.priors = [
            op.join(inia_dir, "inia19-prob_1.nii"),
            op.join(inia_dir, "inia19-prob_2.nii"),
            op.join(inia_dir, "inia19-prob_0.nii")
        ]
    
    main(
        data_path=args.data,
        main_path=args.out,
        subjects=args.subjects,
        sessions=args.ses,
        acquisitions=args.acq,
        template=args.template,
        priors=args.priors if len(args.priors) > 1 else args.priors[0]
    )
