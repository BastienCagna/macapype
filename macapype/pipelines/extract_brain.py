"""
    TODO :p

"""
import os.path as op

import nipype.interfaces.utility as niu
import nipype.pipeline.engine as pe

import nipype.interfaces.fsl as fsl
import nipype.interfaces.afni as afni
import nipype.interfaces.spm as spm

from ..nodes.extract_brain import apply_atlasBREX


def create_brain_extraction_pipe(name="brain_extraction_pipe"):

    # creating pipeline
    brain_extraction_pipe = pe.Workflow(name=name)

    # creating inputnode
    inputnode = pe.Node(
        niu.IdentityInterface(fields=['restore_T1', 'restore_T2']),
        name='inputnode')

    # atlas_brex
    atlas_brex = pe.Node(niu.Function(input_names=['t1_restored_file'],
                                      output_names=['brain_file'],
                                      function=apply_atlasBREX),
                         name='atlas_brex')

    brain_extraction_pipe.connect(inputnode, "restore_T1",
                                  atlas_brex, 't1_restored_file')

    # mask_brex
    mask_brex = pe.Node(fsl.UnaryMaths(), name='mask_brex')
    mask_brex.inputs.operation = 'bin'

    brain_extraction_pipe.connect(atlas_brex, 'brain_file',
                                  mask_brex, 'in_file')

    # smooth_mask
    smooth_mask = pe.Node(fsl.UnaryMaths(), name='smooth_mask')
    smooth_mask.inputs.operation = "bin"
    smooth_mask.inputs.args = "-s 1 -thr 0.5 -bin"

    brain_extraction_pipe.connect(mask_brex, 'out_file',
                                  smooth_mask, 'in_file')

    # mult_T1
    mult_T1 = pe.Node(afni.Calc(), name='mult_T1')
    mult_T1.inputs.expr = "a*b"
    mult_T1.inputs.outputtype = 'NIFTI_GZ'

    brain_extraction_pipe.connect(inputnode, "restore_T1",
                                  mult_T1, 'in_file_a')
    brain_extraction_pipe.connect(smooth_mask, 'out_file',
                                  mult_T1, 'in_file_b')

    # mult_T2
    mult_T2 = pe.Node(afni.Calc(), name='mult_T2')
    mult_T2.inputs.expr = "a*b"
    mult_T2.inputs.outputtype = 'NIFTI_GZ'

    brain_extraction_pipe.connect(inputnode, 'restore_T1',
                                  mult_T2, 'in_file_a')
    brain_extraction_pipe.connect(smooth_mask, 'out_file',
                                  mult_T2, 'in_file_b')
    return brain_extraction_pipe


def create_old_segment_extraction_pipe(name="old_segment_exctraction_pipe"):
    """ Extract brain using tissues masks outputed by SPM's old_segment function

    Inputs
    ---------

    Outputs
    --------
    
    """

    # creating pipeline
    be_pipe = pe.Workflow(name=name)

    # creating inputnode
    inputnode = pe.Node(
        niu.IdentityInterface(fields=['T1'], ['seg_priors']),
        name='inputnode'
    )
    
    # Segment in to 6 tissues
    segment = pe.Node(spm.Segment(), name="old_segment")
    segment.inputs.gm_output_type = [False,False,True]
    segment.inputs.wm_output_type = [False,False,True]
    segment.inputs.csf_output_type = [False,False,True]
    be_pipe.connect(inputnode, 'T1', segment, 'data')
    be_pipe.connect(inputnode, 'seg_priors', segment, 'tissue_prob_maps')

    # Threshold GM, WM and CSF
    thd_nodes = {}
    for tissue in ['gm', 'wm', 'csf']:
        tmp_node = pe.Node(fsl.Threshold(), name="threshold_" + tissue)
        tmp_node.inputs.thresh = 0.05
        be_pipe.connect(
            segment, 'native_' + tissue + '_image', 
            tmp_node, 'in_file'
        )
        thd_nodes[tissue] = tmp_node

    # Compute union of the 3 tissues
    # Done with 2 fslmaths as it seems to hard to do it 
    wmgm_union = pe.Node(fsl.BinaryMaths(), name="wmgm_union")
    wmgm_union.inputs.operation = "add"
    be_pipe.connect(thd_nodes['gm'], 'out_file', wmgm_union, 'in_file')
    be_pipe.connect(thd_nodes['wm'], 'out_file', wmgm_union, 'operand_file')

    tissues_union = pe.Node(fsl.BinaryMaths(), name="wmgm_union")
    tissues_union.inputs.operation = "add"
    be_pipe.connect(wmgm_union, 'out_file', tissues_union, 'in_file')
    be_pipe.connect(thd_nodes['csf'], 'out_file', tissues_union, 'operand_file')
    
    # Opening
    opening_shape = "sphere"
    opening_size = 2
    dilate_mask = pe.Node(fsl.BinaryMaths(), name="dilate_mask")
    # Arbitrary operation
    dilate_mask.inputs.operation = "mean"
    dilate_mask.inputs.kernel_shape = opening_shape
    dilate_mask.inputs.kernel_size = opening_size
    be_pipe.connect(tissues_union, 'out_file', dilate_mask, 'in_file')

    
    erode_mask = pe.Node(fsl.BinaryMaths(), name="erode_mask")
    erode_mask.inputs.kernel_shape = opening_shape
    erode_mask.inputs.kernel_size = opening_size
    be_pipe.connect(dilate_mask, 'out_file', erode_mask, 'in_file')

    # TODO: Add hole filling

    return be_pipe

