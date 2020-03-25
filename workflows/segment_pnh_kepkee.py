#!/usr/bin/env python3
"""
    PNH anatomical segmentation pipeline given by Kepkee Loh wrapped in
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
    python segment_pnh_regis.py -data [PATH_TO_BIDS] -out ../tests/ -subjects Elouk

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
        - ANTS
        - AFNI
        - FSL

"""

# Authors : David Meunier (david.meunier@univ-amu.fr)
#           Bastien Cagna (bastien.cagna@univ-amu.fr)
#           Kepkee Loh (kepkee.loh@univ-amu.fr)
#           Julien Sein (julien.sein@univ-amu.fr)

import os
import os.path as op

import argparse


import nipype

import nipype.interfaces.io as nio
import nipype.interfaces.utility as niu
import nipype.pipeline.engine as pe

import nipype.interfaces.fsl as fsl
fsl.FSLCommand.set_default_output_type('NIFTI_GZ')

#import pybids
from bids.layout import BIDSLayout
from nipype.interfaces.io import BIDSDataGrabber

from macapype.pipelines.preproc import (create_average_align_crop_pipe,
                                        create_average_align_pipe)

from macapype.pipelines.denoise import create_denoised_pipe

from macapype.pipelines.correct_bias import create_correct_bias_pipe
from macapype.pipelines.extract_brain import create_brain_extraction_pipe
from macapype.pipelines.segment import create_full_segment_pipe

from macapype.utils.misc import show_files
from macapype.utils.utils_tests import load_test_data

nmt_dir = load_test_data('NMT_v1.2')
atlasbrex_dir = load_test_data('AtlasBREX')

def create_infosource(subject_ids):
    infosource = pe.Node(interface=niu.IdentityInterface(fields=['subject_id']),name="infosource")
    infosource.iterables = [('subject_id', subject_ids)]

    return infosource


def create_datasource(data_dir, ses):
    datasource = pe.Node(
        interface=nio.DataGrabber(infields=['subject_id'], outfields=['T1', 'T1cropbox', 'T2']),
        name='datasource'
    )
    datasource.inputs.base_directory = data_dir
    #datasource.inputs.template = '%s/sub-%s_ses-01_%s.nii'
    datasource.inputs.template = 'sub-%s/{}/sub-%s_{}_%s.%s'.format(ses, ses)
    datasource.inputs.template_args = dict(
        T1=[['subject_id','subject_id', "mp2rageT1w", "nii"]],
        T1cropbox=[['subject_id', 'subject_id', "mp2rageT1wCropped", "cropbox"]],
        T2=[['subject_id', 'subject_id', "T2w", "nii"]],
    )
    datasource.inputs.sort_filelist = True

    return datasource


def create_datasource_ucdavis_crop_req(data_dir, ses):
    datasource = pe.Node(
        interface=nio.DataGrabber(infields=['subject_id'], outfields=['T1', 'T1cropbox', 'T2']),
        name='datasource'
    )
    datasource.inputs.base_directory = data_dir
    #datasource.inputs.template = '%s/sub-%s_ses-01_%s.nii'
    datasource.inputs.template = 'sub-%s/{}/anat/sub-%s_*run-*_%s.%s'.format(ses)
    datasource.inputs.template_args = dict(
        T1=[['subject_id','subject_id', "T1w", "nii.gz"]],
        T1cropbox=[['subject_id', 'subject_id', "T1wCropped", "cropbox"]],
        T2=[['subject_id', 'subject_id', "T2w", "nii.gz"]],
    )
    datasource.inputs.sort_filelist = True

    return datasource


def create_datasource_ucdavis(data_dir, ses):
    datasource = pe.Node(
        interface=nio.DataGrabber(infields=['subject_id'], outfields=['T1', 'T2']),
        name='datasource'
    )
    datasource.inputs.base_directory = data_dir
    #datasource.inputs.template = '%s/sub-%s_ses-01_%s.nii'
    datasource.inputs.template = 'sub-%s/{}/anat/sub-%s_*run-*_%s.%s'.format(ses)
    datasource.inputs.template_args = dict(
        T1=[['subject_id','subject_id', "T1w", "nii.gz"]],
        T2=[['subject_id', 'subject_id', "T2w", "nii.gz"]],
    )
    datasource.inputs.sort_filelist = True

    return datasource


def create_bids_datasource(data_dir):

    bids_datasource = pe.Node(
        interface=nio.BIDSDataGrabber(),
        name='bids_datasource'
    )

    bids_datasource.inputs.base_dir = data_dir
    bids_datasource.inputs.output_query = {
        'T1': {"datatype": "anat",
               "suffix": "T1w",
               "extensions": ["nii", ".nii.gz"]},
        'T2': {"datatype": "anat",
               "suffix": "T2w",
               "extensions": ["nii", ".nii.gz"]}}

    layout = BIDSLayout(data_dir)
    print (layout)
    print(layout.get_subjects())
    print(layout.get_sessions())


    iterables = []

    assert len(layout.get_subjects()) or len(layout.get_sessions()), \
        "Error, BIDS dir should have at least one subject and one session"

    if len(layout.get_subjects()) == 1:
         bids_datasource.inputs.subject = layout.get_subjects()[0]
    else:
        iterables.append(('subject', layout.get_subjects()))

    if len(layout.get_sessions()) == 1:
         bids_datasource.inputs.session = layout.get_sessions()[0]
    else:
        iterables.append(('session', layout.get_sessions()))

    if len(iterables):
        bids_datasource.iterables = iterables

    return bids_datasource

###############################################################################

def create_segment_pnh_subpipes(crop_req, name= "segment_pnh_subpipes",
                                sigma = 2):

    # creating pipeline
    seg_pipe = pe.Workflow(name=name)

    """
    new version (as it is now)
    - preproc (avg and align, crop is optional)
    - correct_bias
    - denoise
    - extract_brain
    - segment
    """
    print("crop_req:", crop_req)

    if crop_req:

        inputnode = pe.Node(
            niu.IdentityInterface(fields=['T1', 'T1cropbox', 'T2']),
            name='inputnode')

        #### preprocessing (avg and align)
        preproc_pipe = create_average_align_crop_pipe()

        seg_pipe.connect(inputnode,'T1',preproc_pipe,'inputnode.T1')
        seg_pipe.connect(inputnode,'T2',preproc_pipe,'inputnode.T2')

        seg_pipe.connect(inputnode,'T1cropbox',preproc_pipe,'inputnode.T1cropbox')
        seg_pipe.connect(inputnode,'T1cropbox',preproc_pipe,'inputnode.T2cropbox')

        #### Correct_bias_T1_T2
        correct_bias_pipe = create_correct_bias_pipe(sigma = sigma)

        seg_pipe.connect(preproc_pipe, "crop_bb_T1.roi_file",correct_bias_pipe,'inputnode.preproc_T1')
        seg_pipe.connect(preproc_pipe, "crop_bb_T2.roi_file", correct_bias_pipe,'inputnode.preproc_T2')

        ##### otherwise using nibabel node
        #from nodes.segment_pnh_nodes import correct_bias_T1_T2

        #correct_bias = pe.Node(interface = niu.Function(
            #input_names=["preproc_T1_file","preproc_T2_file", "sigma"],
            #output_names =  ["thresh_lower_file", "norm_mult_file", "bias_file", "smooth_bias_file", "restore_T1_file", "restore_T2_file"],
            #function = correct_bias_T1_T2),
            #name = "correct_bias")

        #correct_bias.inputs.sigma = sigma*2

        #seg_pipe.connect(preproc_pipe, 'crop_bb_T1.roi_file',correct_bias,'preproc_T1_file')
        #seg_pipe.connect(preproc_pipe, 'crop_bb_T2.roi_file',correct_bias,'preproc_T2_file')

    else:
        inputnode = pe.Node(
            niu.IdentityInterface(fields=['T1', 'T2']),
            name='inputnode')

        #### preprocessing (avg and align)
        preproc_pipe = create_average_align_pipe()

        seg_pipe.connect(inputnode,'T1',preproc_pipe,'inputnode.T1')
        seg_pipe.connect(inputnode,'T2',preproc_pipe,'inputnode.T2')

        #### Correct_bias_T1_T2
        correct_bias_pipe = create_correct_bias_pipe(sigma = sigma)

        seg_pipe.connect(preproc_pipe, "av_T1.avg_img", correct_bias_pipe,'inputnode.preproc_T1')
        seg_pipe.connect(preproc_pipe, "align_T2_on_T1.out_file", correct_bias_pipe,'inputnode.preproc_T2')


    #### denoising
    denoise_pipe = create_denoised_pipe()

    seg_pipe.connect(correct_bias_pipe, "restore_T1.out_file", denoise_pipe,'inputnode.preproc_T1')
    seg_pipe.connect(correct_bias_pipe, "restore_T2.out_file",denoise_pipe,'inputnode.preproc_T2')

    ##### brain extraction
    brain_extraction_pipe = create_brain_extraction_pipe(
        atlasbrex_dir=atlasbrex_dir, nmt_dir=nmt_dir, name = "devel_atlas_brex")

    seg_pipe.connect(denoise_pipe,'denoise_T1.output_image',brain_extraction_pipe,"inputnode.restore_T1")
    seg_pipe.connect(denoise_pipe,'denoise_T2.output_image',brain_extraction_pipe,"inputnode.restore_T2")

    ################### full_segment (restarting from the avg_align files,)
    brain_segment_pipe = create_full_segment_pipe(sigma=sigma, nmt_dir = nmt_dir,
        name="segment_devel_NMT_sub_align")

    #cropped False
    seg_pipe.connect(preproc_pipe, "av_T1.avg_img" ,brain_segment_pipe,'inputnode.preproc_T1')
    seg_pipe.connect(preproc_pipe, "align_T2_on_T1.out_file", brain_segment_pipe,'inputnode.preproc_T2')

    # Cropped True
    #seg_pipe.connect(preproc_pipe, "crop_bb_T1.roi_file",brain_segment_pipe,'inputnode.preproc_T1')
    #seg_pipe.connect(preproc_pipe, "crop_bb_T2.roi_file", brain_segment_pipe,'inputnode.preproc_T2')

    seg_pipe.connect(brain_extraction_pipe,"smooth_mask.out_file",brain_segment_pipe,"inputnode.brain_mask")

    return seg_pipe

def create_main_workflow(data_dir, process_dir, subject_ids, ses, crop_req):

    main_workflow = pe.Workflow(name= "test_pipeline_david_no_crop_req")
    main_workflow.base_dir = process_dir


    print ("crop_req:", crop_req)


    #if subject_ids is None or sess is None and cropped is not None:
        #print('adding BIDS data source')
        #datasource = create_bids_datasource(data_dir)
        #print(datasource.outputs)

    #else:

    print('adding infosource and datasource')

    # Infosource
    infosource = create_infosource(subject_ids)

    # Data source
    if crop_req is True:
        print ("Datasource crop_req")
        datasource = create_datasource_ucdavis_crop_req(data_dir, ses)
    else:
        print ("Datasource already cropped")
        datasource = create_datasource_ucdavis(data_dir, ses)

    # connect
    main_workflow.connect(infosource, 'subject_id', datasource, 'subject_id')

    ############################################## Preprocessing ################################
    ##### segment_pnh

    print('segment_pnh')

    segment_pnh = create_segment_pnh_subpipes(crop_req = crop_req)

    main_workflow.connect(datasource,'T1',segment_pnh,'inputnode.T1')
    main_workflow.connect(datasource,'T2',segment_pnh,'inputnode.T2')

    if crop_req is True:
        main_workflow.connect(datasource,'T1cropbox',segment_pnh,'inputnode.T1cropbox')

    return main_workflow

if __name__ == '__main__':

    # Command line parser
    parser = argparse.ArgumentParser(
        description="PNH segmentation pipeline from Kepkee Loh / Julien Sein")

    parser.add_argument("-data", dest="data", type=str, required=True,
                        help="Directory containing MRI data (BIDS)")
    parser.add_argument("-out", dest="out", type=str, #nargs='+',
                        help="Output dir", required=True)
    parser.add_argument("-ses", dest="ses", type=str,
                        help="Session", required=False)
    parser.add_argument("-subjects", dest="subjects", type=str, nargs='+',
                        help="Subjects' ID", required=False)
    parser.add_argument("-crop_req", dest="crop_req", default = False,
                        action='store_true', help="Crop required")

    args = parser.parse_args()

    # main_workflow
    print("Initialising the pipeline...")
    wf = create_main_workflow(
        data_dir=args.data,
        process_dir=args.out,
        subject_ids=args.subjects,
        ses=args.ses,
        crop_req=args.crop_req
    )
    wf.write_graph(graph2use="colored")
    wf.config['execution'] = {'remove_unnecessary_outputs': 'false'}
    print('The PNH segmentation pipeline is ready')

    print("Start to process")
    wf.run()
    # wf.run(plugin='MultiProc', plugin_args={'n_procs' : 2})