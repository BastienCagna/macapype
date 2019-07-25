import nipype

import nipype.interfaces.io as nio
import nipype.interfaces.utility as niu
import nipype.pipeline.engine as pe

import nipype.interfaces.fsl as fsl
fsl.FSLCommand.set_default_output_type('NIFTI_GZ')

import os

from macapype.pipelines.segment import create_full_segment_pipe

from macapype.utils.misc import show_files

#data_path = "/hpc/meca/data/Macaques/Macaque_hiphop/"
#main_path = "/hpc/crise/meunier.d/Data/"
main_path = "/hpc/meca/users/loh.k/test_pipeline"
site = "sbri"
subject_ids = ['032311']

prev_pipe_path = os.path.join(main_path, "test_pipeline_kepkee_by_kepkee",
                              "segment_pnh_subpipes")

def create_infosource():
    infosource = pe.Node(interface=niu.IdentityInterface(fields=['subject_id']),name="infosource")
    infosource.iterables = [('subject_id', subject_ids)]

    return infosource

def create_datasource():
   datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],outfields=['T1','T2', 'mask']),name = 'datasource')
   datasource.inputs.base_directory = prev_pipe_path
   datasource.inputs.template = '%s/%s%s/%s%s/%ssub-%s_ses-001_run-*_%s.nii.gz'
   datasource.inputs.template_args = dict(
       T1=[["average_align_pipe","_subject_id_", 'subject_id',"","av_T1","",
            'subject_id',"T1w_flirt"]],
       T2=[["average_align_pipe","_subject_id_", 'subject_id',"","align_T2_on_T1",
            "",'subject_id',"T2w_flirt"]],
       mask=[["brain_extraction_pipe", "_subject_id_",
              'subject_id', "","smooth_mask", "avg_", 'subject_id',
              #'T1w_maths_aonlm_denoised_roi_brain_bin_bin_nice']],
              'T1w_maths_aonlm_denoised_roi_brain_bin_bin']],
       )
   datasource.inputs.sort_filelist = True

   return datasource

#def create_datasource():
   #datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],outfields=['T1','T2']),name = 'datasource')
   #datasource.inputs.base_directory = "/hpc/meca/"
   #datasource.inputs.template = '%s/%s/%s%s/%s/%s/sub-%s_ses-001_run-*_%s.nii.gz'
   #datasource.inputs.template_args = dict(
       #T1=[["data/Macaques/Macaque_hiphop/",site,"sub-",'subject_id',"ses-001","anat",'subject_id',"T1w"]],
       #T2=[["data/Macaques/Macaque_hiphop/",site,"sub-",'subject_id',"ses-001","anat",'subject_id',"T2w"]],
       #mask=[["users/loh.k/test_pipeline","test_pipeline_kepkee_by_kepkee","segment_pnh_subpipes","brain_extraction_pipe","smooth_mask"
              #site,'subject_id',"mask",'subject_id',"mask"]],
       #)
   #datasource.inputs.sort_filelist = True

   #return datasource

###############################################################################

def create_main_workflow():

    main_workflow = pe.Workflow(name= "test_pipeline_full_segment")
    main_workflow.base_dir = main_path

    ## Infosource
    infosource = create_infosource()

    ## Data source
    datasource = create_datasource()

    ## connect
    main_workflow.connect(infosource, 'subject_id', datasource, 'subject_id')

    ############################################## Preprocessing ################################
    ##### segment_pnh

    print('full_segment_pnh')

    ################### full_segment
    brain_segment_pipe = create_full_segment_pipe(crop_list = [(88, 144), (14, 180), (27, 103)], sigma = 2)

    main_workflow.connect(datasource, "T1",
                          brain_segment_pipe, 'inputnode.preproc_T1')
    main_workflow.connect(datasource, "T2",
                          brain_segment_pipe, 'inputnode.preproc_T2')

    main_workflow.connect(datasource, "mask",brain_segment_pipe,"inputnode.brain_mask")

    return main_workflow


if __name__ =='__main__':

    ### main_workflow
    wf = create_main_workflow()
    wf.write_graph(graph2use = "colored")
    wf.config['execution'] = {'remove_unnecessary_outputs':'false'}

    wf.run()
    #wf.run(plugin='MultiProc', plugin_args={'n_procs' : 2})
