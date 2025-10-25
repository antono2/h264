module h264

// Minimal H264 video parser
// Read the H264 specification: Rec. ITU-T H.264 (08/2021)
// The following library was used as reference: https://github.com/aizvorski/h264bitstream/
// 
// Use this to parse metadata information about a H264 video before passing it to a decoder.
// You can use the following functions to access H264 information:
//
// Create byte array to consume binary data:
// mut nal := h264.NetworkAbstractionLayerHeader{}
// mut nal_header_bs := h264.Bitstream{}
// nal_header_bs.init(data_sps.vbytes(1))
//
// Then you can read a NAL header after you detected a new NAL unit (starts with NalStartCode bytes):
// nal.read_nal_header(mut nal_header_bs)
//
// After the NAL unit type was determined, you can use the following functions to read PPS, SPS and SliceHeaders (depending on NAL unit type):
// pub fn (mut pps PictureParameterSet) read_pps(mut b Bitstream)
// pub fn (mut sps SequenceParameterSet) read_sps(mut b Bitstream)
// pub fn (mut sh SliceHeader) read_slice_header(nal &NetworkAbstractionLayerHeader, pps_array []PictureParameterSet, sps_array []SequenceParameterSet, mut b Bitstream)
//


//pub const nal_start_code = [u8(0), 0, 1]
pub struct NalStartCode {
pub:
  value []u8 = [u8(0), 0, 1]
}

pub struct SequenceParameterSet {
pub mut:
  profile_idc u32
  constraint_set0_flag u32
  constraint_set1_flag u32
  constraint_set2_flag u32
  constraint_set3_flag u32
  constraint_set4_flag u32
  constraint_set5_flag u32
  reserved_zero_2bits i32
  level_idc u32
  seq_parameter_set_id u32
  chroma_format_idc u32
  separate_colour_plane_flag u32
  bit_depth_luma_minus8 u32
  bit_depth_chroma_minus8 u32
  qpprime_y_zero_transform_bypass_flag u32
  seq_scaling_matrix_present_flag u32
  seq_scaling_list_present_flag [8]u32
  scaling_list_4x4 [6][16]i32
  use_default_scaling_matrix_4x4_flag [6]u32
  scaling_list_8x8 [2][64]i32
  use_default_scaling_matrix_8x8_flag [2]u32
  log2_max_frame_num_minus4 u32
  pic_order_cnt_type u32
  log2_max_pic_order_cnt_lsb_minus4 u32
  delta_pic_order_always_zero_flag u32
  offset_for_non_ref_pic i32
  offset_for_top_to_bottom_field i32
  num_ref_frames_in_pic_order_cnt_cycle u32
  offset_for_ref_frame [256]i32
  num_ref_frames u32
  gaps_in_frame_num_value_allowed_flag u32
  pic_width_in_mbs_minus1 u32
  pic_height_in_map_units_minus1 u32
  frame_mbs_only_flag u32
  mb_adaptive_frame_field_flag u32
  direct_8x8_inference_flag u32
  frame_cropping_flag u32
  frame_crop_left_offset u32
  frame_crop_right_offset u32
  frame_crop_top_offset u32
  frame_crop_bottom_offset u32
  vui_parameters_present_flag u32
  vui SequenceParameterSetVideoUsabilityInformation
  hrd SequenceParameterSetHypotheticalReferenceDecoder
}

pub struct SequenceParameterSetVideoUsabilityInformation {
pub mut:
  aspect_ratio_info_present_flag u32
  aspect_ratio_idc u32
  sar_width u32
  sar_height u32
  overscan_info_present_flag u32
  overscan_appropriate_flag u32
  video_signal_type_present_flag u32
  video_format u32
  video_full_range_flag u32
  color_description_present_flag u32
  colour_primaries u32
  transfer_characteristics u32
  matrix_coefficients u32
  chroma_loc_info_present_flag u32
  chroma_sample_loc_type_top_field u32
  chroma_sample_loc_type_bottom_field u32
  timing_info_present_flag u32
  num_units_in_tick u32
  time_scale u32
  fixed_frame_rate_flag u32
  nal_hrd_parameters_present_flag u32
  vcl_hrd_parameters_present_flag u32
  low_delay_hrd_flag u32
  pic_struct_present_flag u32
  bitstream_restriction_flag u32
  motion_vectors_over_pic_boundaries_flag u32
  max_bytes_per_pic_denom u32
  max_bits_per_mb_denom u32
  log2_max_mv_length_horizontal u32
  log2_max_mv_length_vertical u32
  num_reorder_frames u32
  max_dec_frame_buffering u32
}

pub struct SequenceParameterSetHypotheticalReferenceDecoder {
pub mut:
  cpb_cnt_minus1 u32
  bit_rate_scale u32
  cpb_size_scale u32
  bit_rate_value_minus1 [32]u32
  cpb_size_value_minus1 [32]u32
  cbr_flag [32]u32
  initial_cpb_removal_delay_length_minus1 u32
  cpb_removal_delay_length_minus1 u32
  dpb_output_delay_length_minus1 u32
  time_offset_length u32
}

pub  struct PictureParameterSet {
pub mut:
  pic_parameter_set_id u32
  seq_parameter_set_id u32
  entropy_coding_mode_flag u32
  pic_order_present_flag u32
  num_slice_groups_minus1 u32
  slice_group_map_type u32
  run_length_minus1 [8]u32
  top_left [8]u32
  bottom_right [8]u32
  slice_group_change_direction_flag u32
  slice_group_change_rate_minus1 u32
  pic_size_in_map_units_minus1 u32
  slice_group_id [256]u32
  num_ref_idx_l0_active_minus1 u32
  num_ref_idx_l1_active_minus1 u32
  weighted_pred_flag u32
  weighted_bipred_idc u32
  pic_init_qp_minus26 i32
  pic_init_qs_minus26 i32
  chroma_qp_index_offset i32
  deblocking_filter_control_present_flag u32
  constrained_intra_pred_flag u32
  redundant_pic_cnt_present_flag u32
  is_more_rbsp_data_present i32
  transform_8x8_mode_flag u32
  pic_scaling_matrix_present_flag u32
  pic_scaling_list_present_flag [8]u32
  scaling_list_4x4 [6][16]i32
  use_default_scaling_matrix_4x4_flag [6]u32
  scaling_list_8x8 [2][64]i32
  use_default_scaling_matrix_8x8_flag [2]u32
  second_chroma_qp_index_offset i32
}

pub enum SH_SLICE_TYPE {
  p = 0
  b = 1
  i = 2
  sp = 3
  si = 4
  // The ONLY slice types indicating that all other slices in that picture are of the same type
  p_only = 5
  b_only = 6
  i_only = 7
  sp_only = 8
  si_only = 9
}

pub struct SliceHeader {
pub mut:
  first_mb_in_slice u32
  slice_type u32
  pic_parameter_set_id u32
  frame_num u32
  field_pic_flag u32
  bottom_field_flag u32
  idr_pic_id u32
  pic_order_cnt_lsb u32
  delta_pic_order_cnt_bottom i32
  delta_pic_order_cnt [2]i32
  redundant_pic_cnt u32
  direct_spatial_mv_pred_flag u32
  num_ref_idx_active_override_flag u32
  num_ref_idx_l0_active_minus1 u32
  num_ref_idx_l1_active_minus1 u32
  cabac_init_idc u32
  slice_qp_delta i32
  sp_for_switch_flag u32
  slice_qs_delta i32
  disable_deblocking_filter_idc u32
  slice_alpha_c0_offset_div2 i32
  slice_beta_offset_div2 i32
  slice_group_change_cycle u32
  pwt SliceHeaderPredictiveWeightTable
  rplr SliceHeaderReferencePictureListReorder
  drpm SliceHeaderDecodedReferencePictureMarking
}

pub struct SliceHeaderPredictiveWeightTable {
pub mut:
  luma_log2_weight_denom u32
  chroma_log2_weight_denom u32
  luma_weight_l0_flag [64]u32
  luma_weight_l0 [64]i32
  luma_offset_l0 [64]i32
  chroma_weight_l0_flag [64]u32
  chroma_weight_l0 [64][2]i32
  chroma_offset_l0 [64][2]i32
  luma_weight_l1_flag [64]u32
  luma_weight_l1 [64]i32
  luma_offset_l1 [64]i32
  chroma_weight_l1_flag [64]u32
  chroma_weight_l1 [64][2]i32
  chroma_offset_l1 [64][2]i32
}

pub struct SliceHeaderReferencePictureListReorder {
pub mut:
  ref_pic_list_reordering_flag_l0 u32
  reorder_l0 SliceHeaderReferencePictureListReorderL0
  ref_pic_list_reordering_flag_l1 u32
  reorder_l1 SliceHeaderReferencePictureListReorderL1
}

pub struct SliceHeaderReferencePictureListReorderL0 {
pub mut:
  reordering_of_pic_nums_idc [64]u32
  abs_diff_pic_num_minus1 [64]u32
  long_term_pic_num [64]u32
}

pub struct SliceHeaderReferencePictureListReorderL1 {
pub mut:
  reordering_of_pic_nums_idc [64]u32
  abs_diff_pic_num_minus1 [64]u32
  long_term_pic_num [64]u32
}

pub struct SliceHeaderDecodedReferencePictureMarking {
pub mut:
  no_output_of_prior_pics_flag u32
  long_term_reference_flag u32
  adaptive_ref_pic_marking_mode_flag u32
  memory_management_control_operation [64]u32
  difference_of_pic_nums_minus1 [64]u32
  long_term_pic_num [64]u32
  long_term_frame_idx [64]u32
  max_long_term_frame_idx_plus1 [64]u32
}

pub enum NAL_REF_IDC {
  priority_highest = 3
  priority_high = 2
  priority_low = 1
  priority_disposable = 0
}

pub enum NAL_UNIT_TYPE {
  unspecified = 0	// Unspecified
  coded_slice_non_idr = 1	// Coded slice of a non-IDR picture
  coded_slice_data_partition_a = 2	// Coded slice data partition A
  coded_slice_data_partition_b = 3	// Coded slice data partition B
  coded_slice_data_partition_c = 4	// Coded slice data partition C
  coded_slice_idr = 5	// Coded slice of an IDR picture
  sei = 6	// Supplemental enhancement information (SEI)
  sps = 7	// Sequence parameter set
  pps = 8	// Picture parameter set
  aud = 9	// Access unit delimiter
  end_of_sequence = 10	// End of sequence
  end_of_stream = 11	// End of stream
  filler = 12	// Filler data
  sps_ext = 13	// Sequence parameter set extension
  coded_slice_aux = 19	// Coded slice of an auxiliary coded picture without partitioning
}

pub struct NetworkAbstractionLayerHeader {
pub mut:
  idc NAL_REF_IDC
  type NAL_UNIT_TYPE
}

pub struct Bitstream {
pub mut:
  start int
  p int
  end int
  b []u8
  bits_left u8
}

// Use byte array to not deal with pointers too much.
pub fn (mut b Bitstream) init(buf []u8) {
  b.b = buf
  b.p = 0
  b.end = buf.len
  b.bits_left = 8
}

pub fn (b Bitstream) byte_aligned() bool {
  return b.bits_left == 8
} 

pub fn (b Bitstream) eof() bool {
  if b.p >= b.end {
    return true
  }
  return false
}

// Read next bit and keep track of bits left.
// Then jump to the next byte after 0 bits left
pub fn (mut b Bitstream) u1() u32 {
  mut r := u32(0)
  b.bits_left--
  if !b.eof() {
    r = (b.b[b.p] >> b.bits_left) & 0x01
  }
  if b.bits_left == 0 {
    b.p++
    b.bits_left = 8
  }
  return r
}

// Read unsigned int, n bits long
pub fn (mut b Bitstream) u(n u32) u32 {
  mut r := u32(0)
  for i in 0..n {
    r |= b.u1() << (n - i - 1)
  }
  return r
}

// Read unsigned Exp-Golomb number
// Count the bits until 1, that's the number of following bits to read into r
pub fn (mut b Bitstream) ue() u32 {
  mut i := u32(0)
  
  for (b.u1() == 0) && (i < 32) && (!b.eof()) {
    i++
  }
  mut r := b.u(i)
  r += (1 << i) - 1
  return r
}

// Read a signed Exp-Golomb number
pub fn (mut b Bitstream) se() i32 {
  mut r := i32(b.ue())
  if (r & 0x01) != 0 {
    r = (r + 1) / 2
  } else {
    r = -(r / 2)
  }
  return r
}

pub fn (mut nal NetworkAbstractionLayerHeader) read_nal_header(mut b Bitstream) {
  forbidden_zero_bit := b.u1()
  assert forbidden_zero_bit == 0
  nal.idc = unsafe{ NAL_REF_IDC(b.u(2)) }
  nal.type = unsafe{ NAL_UNIT_TYPE(b.u(5)) }
}

pub fn (mut b Bitstream) read_scaling_list(mut scaling_list []i32, size_of_scaling_list i32, mut use_default_scaling_matrix_flag &u32) {
  mut last_scale := i32(8)
  mut next_scale := i32(8)
  mut delta_scale := i32(0)

  for j in 0..size_of_scaling_list {
    if next_scale != 0 {
      delta_scale = b.se()
      next_scale = (last_scale + delta_scale + 256) % 256
      if j == 0 && next_scale == 0 {
        unsafe{*use_default_scaling_matrix_flag = 1}
      } else {
        unsafe{*use_default_scaling_matrix_flag = 0}
      }
    }
    unsafe {
      scaling_list[j] = if next_scale == 0 { last_scale } else { next_scale }
      last_scale = scaling_list[j]
    }
  }
}

pub fn (mut sps SequenceParameterSet) read_hrd_parameters(mut b Bitstream) {
  sps.hrd.cpb_cnt_minus1 = b.ue()
  sps.hrd.bit_rate_scale = b.u(4)
  sps.hrd.cpb_size_scale = b.u(4)
  for sched_sel_idx in 0 .. sps.hrd.cpb_cnt_minus1 {
    sps.hrd.bit_rate_value_minus1[sched_sel_idx] = b.ue()
    sps.hrd.cpb_size_value_minus1[sched_sel_idx] = b.ue()
    sps.hrd.cbr_flag[sched_sel_idx] = b.u1()
  }
  sps.hrd.initial_cpb_removal_delay_length_minus1 = b.u(5)
  sps.hrd.cpb_removal_delay_length_minus1 = b.u(5)
  sps.hrd.dpb_output_delay_length_minus1 = b.u(5)
  sps.hrd.time_offset_length = b.u(5)
}

pub fn (mut b Bitstream) read_rbsp_trailing_bits() {
  // rbsp_stop_one_bit
  b.u1()
  for !b.byte_aligned() {
    // rbsp_alignment_zero_bit
    b.u1()
  }
}

pub fn (mut sps SequenceParameterSet) read_vui_parameters(mut b Bitstream) {
  dump(b)
  sps.vui.aspect_ratio_info_present_flag = b.u1()
  if sps.vui.aspect_ratio_info_present_flag != 0 {
    sps.vui.aspect_ratio_idc = b.u(8)
    // Extended SAR
    if sps.vui.aspect_ratio_idc == 255 {
      sps.vui.sar_width = b.u(16)
      sps.vui.sar_height = b.u(16)
    }
  }
  sps.vui.overscan_info_present_flag = b.u1()
  if sps.vui.overscan_info_present_flag != 0 {
    sps.vui.overscan_appropriate_flag = b.u1()
  }
  sps.vui.video_signal_type_present_flag = b.u1()
  if sps.vui.video_signal_type_present_flag != 0 {
    sps.vui.video_format = b.u(3)
    sps.vui.video_full_range_flag = b.u1()
    sps.vui.color_description_present_flag = b.u1()
    if sps.vui.color_description_present_flag != 0 {
      sps.vui.colour_primaries = b.u(8)
      sps.vui.transfer_characteristics = b.u(8)
      sps.vui.matrix_coefficients = b.u(8)
    }
  }
  sps.vui.chroma_loc_info_present_flag = b.u1()
  if sps.vui.chroma_loc_info_present_flag != 0 {
    sps.vui.chroma_sample_loc_type_top_field = b.ue()
    dump(sps.vui.chroma_sample_loc_type_top_field)
    sps.vui.chroma_sample_loc_type_bottom_field = b.ue()
    dump(sps.vui.chroma_sample_loc_type_bottom_field)
  }
  sps.vui.timing_info_present_flag = b.u1()
  if sps.vui.timing_info_present_flag != 0 {
    sps.vui.num_units_in_tick = b.u(32)
    sps.vui.time_scale = b.u(32)
    sps.vui.fixed_frame_rate_flag = b.u1()
  }
  sps.vui.nal_hrd_parameters_present_flag = b.u1()
  if sps.vui.nal_hrd_parameters_present_flag != 0 {
    sps.read_hrd_parameters(mut b)
  }
  sps.vui.vcl_hrd_parameters_present_flag = b.u1()
  if sps.vui.vcl_hrd_parameters_present_flag != 0 {
    sps.read_hrd_parameters(mut b)
  }
  if sps.vui.nal_hrd_parameters_present_flag != 0 || sps.vui.vcl_hrd_parameters_present_flag != 0 {
    sps.vui.low_delay_hrd_flag = b.u1()
  }
  sps.vui.pic_struct_present_flag = b.u1()
  sps.vui.bitstream_restriction_flag = b.u1()
  if sps.vui.bitstream_restriction_flag != 0 {
    sps.vui.motion_vectors_over_pic_boundaries_flag = b.u1()
    sps.vui.max_bytes_per_pic_denom = b.ue()
    sps.vui.max_bits_per_mb_denom = b.ue()
    sps.vui.log2_max_mv_length_horizontal = b.ue()
    sps.vui.log2_max_mv_length_vertical = b.ue()
    sps.vui.num_reorder_frames = b.ue()
    sps.vui.max_dec_frame_buffering = b.ue()
  }
}


pub fn intlog2(x i32) i32 {
  mut log := i32(0)
  mut xx := x
  if xx < 0 { xx = 0 }
  for (xx >> log) > 0 {
    log++
  }
  if log > 0 && xx == (1 << (log - 1)) {
    log--
  }
  return log
}

/*
pub const multiply_de_bruijn_bit_position = 
[
  u32(0), 9, 1, 10, 13, 21, 2, 29, 11, 14, 16, 18, 22, 25, 3, 30,
  8, 12, 20, 28, 15, 17, 24, 7, 19, 27, 23, 6, 26, 5, 4, 31
]

pub fn intlog2(val i32) i32 {
  mut v := val // Make mutable
  v |= v >> 1
  v |= v >> 2
  v |= v >> 4
  v |= v >> 8
  v |= v >> 16
  return i32(multiply_de_bruijn_bit_position[(v * 0x07C4ACDD) >> 27])
}
*/

pub fn (mut b Bitstream) more_rbsp_data() bool {
  // No more data
  if b.eof() { return false }
  
  // Don't copy the whole stream, just restore the original index
  /*
  mut bs_tmp := Bitstream{
    ...
    b
  }

  // No rbsp_stop_bit yet
  if bs_tmp.u1() == 0 { return true }
  
  for !bs_tmp.eof() {
    // A later bit was 1, it wasn't the rsbp_stop_bit
    if bs_tmp.u1() == 1 { return true }
  }
  */
  bk_index := b.p
  defer { b.p = bk_index }
  
  // No rbsp_stop_bit yet
  if b.u1() == 0 { return true }
  
  for !b.eof() {
    // A later bit was 1, it wasn't the rsbp_stop_bit
    if b.u1() == 1 { return true }
  }
  
  // All following bits were 0, it was the rsbp_stop_bit
  return false
}

pub fn (mut sps SequenceParameterSet) read_sps(mut b Bitstream) {
  sps.profile_idc = b.u(8)
  sps.constraint_set0_flag = b.u1()
  sps.constraint_set1_flag = b.u1()
  sps.constraint_set2_flag = b.u1()
  sps.constraint_set3_flag = b.u1()
  sps.constraint_set4_flag = b.u1()
  sps.constraint_set5_flag = b.u1()
  // reserved_zero_2bits
  b.u(2)
  sps.level_idc = b.u(8)
  sps.seq_parameter_set_id = b.ue()
  dump(sps.seq_parameter_set_id)
  
  if sps.profile_idc == 100
  || sps.profile_idc == 110
  || sps.profile_idc == 122
  || sps.profile_idc == 144 {
    sps.chroma_format_idc = b.ue()
    if sps.chroma_format_idc == 3 {
      sps.separate_colour_plane_flag = b.u1()
    }
    sps.bit_depth_luma_minus8 = b.ue()
    sps.bit_depth_chroma_minus8 = b.ue()
    sps.qpprime_y_zero_transform_bypass_flag = b.u1()
    sps.seq_scaling_matrix_present_flag = b.u1()
    if sps.seq_scaling_matrix_present_flag != 0 {
      for i in 0..8 {
        sps.seq_scaling_list_present_flag[i] = b.u1()
        if sps.seq_scaling_list_present_flag[i] != 0 {
          if i < 6 {
            b.read_scaling_list(mut voidptr(&sps.scaling_list_4x4[i]), 16, mut &sps.use_default_scaling_matrix_4x4_flag[i])
          } else {
            b.read_scaling_list(mut voidptr(sps.scaling_list_8x8[i - 6]), 64, mut &sps.use_default_scaling_matrix_8x8_flag[i - 6])
          }
        }
      }
    }
  } // if sps.profile_idc == 100 ...
  sps.log2_max_frame_num_minus4 = b.ue()
  sps.pic_order_cnt_type = b.ue()
  if sps.pic_order_cnt_type == 0 {
    sps.log2_max_pic_order_cnt_lsb_minus4 = b.ue()
  } else {
    if sps.pic_order_cnt_type == 1 {
      sps.delta_pic_order_always_zero_flag = b.u1()
      sps.offset_for_non_ref_pic = b.se()
      sps.offset_for_top_to_bottom_field = b.se()
      sps.num_ref_frames_in_pic_order_cnt_cycle = b.ue()
      for i in 0..sps.num_ref_frames_in_pic_order_cnt_cycle {
        sps.offset_for_ref_frame[i] = b.se()
      }
    }
  }
  sps.num_ref_frames = b.ue()
  sps.gaps_in_frame_num_value_allowed_flag = b.u1()
  sps.pic_width_in_mbs_minus1 = b.ue()
  sps.pic_height_in_map_units_minus1 = b.ue()
  sps.frame_mbs_only_flag = b.u1()
  if sps.frame_mbs_only_flag == 0 {
    sps.mb_adaptive_frame_field_flag = b.u1()
  }
  sps.direct_8x8_inference_flag = b.u1()
  sps.frame_cropping_flag = b.u1()
  if sps.frame_cropping_flag != 0 {
    sps.frame_crop_left_offset = b.ue()
    sps.frame_crop_right_offset = b.ue()
    sps.frame_crop_top_offset = b.ue()
    sps.frame_crop_bottom_offset = b.ue()
  }
  sps.vui_parameters_present_flag = b.u1()
  if sps.vui_parameters_present_flag != 0 {
    sps.read_vui_parameters(mut b)
  }
  b.read_rbsp_trailing_bits()
}

pub fn (mut pps PictureParameterSet) read_pps(mut b Bitstream) {
  pps.pic_parameter_set_id = b.ue()
  dump(pps.pic_parameter_set_id)
  pps.seq_parameter_set_id = b.ue()
  dump(pps.seq_parameter_set_id)
  pps.entropy_coding_mode_flag = b.u1()
  pps.pic_order_present_flag = b.u1()
  pps.num_slice_groups_minus1 = b.ue()

  if pps.num_slice_groups_minus1 > 0 {
    pps.slice_group_map_type = b.ue()
    if pps.slice_group_map_type == 0 {
      for i_group in 0..pps.num_slice_groups_minus1 {
        pps.run_length_minus1[i_group] = b.ue()
      }
    } else if pps.slice_group_map_type == 2 {
      for i_group in 0..pps.num_slice_groups_minus1 {
        pps.top_left[i_group] = b.ue()
        pps.bottom_right[i_group] = b.ue()
      }
    } else if pps.slice_group_map_type == 3
           || pps.slice_group_map_type == 4
           || pps.slice_group_map_type == 5 {
      pps.slice_group_change_direction_flag = b.u1()
      pps.slice_group_change_rate_minus1 = b.ue()
    } else if pps.slice_group_map_type == 6 {
      pps.pic_size_in_map_units_minus1 = b.ue()
      for i in 0..pps.pic_size_in_map_units_minus1 {
        v := intlog2(i32(pps.num_slice_groups_minus1) + 1)
        pps.slice_group_id[i] = b.u(u32(v))
      }
    }
  } // if pps.num_slice_groups_minus1 > 0
  pps.num_ref_idx_l0_active_minus1 = b.ue()
  dump(pps.num_ref_idx_l0_active_minus1)
  pps.num_ref_idx_l1_active_minus1 = b.ue()
  pps.weighted_pred_flag = b.u1()
  pps.weighted_bipred_idc = b.u(2)
  pps.pic_init_qp_minus26 = b.se()
  pps.pic_init_qs_minus26 = b.se()
  pps.chroma_qp_index_offset = b.se()
  pps.deblocking_filter_control_present_flag = b.u1()
  pps.constrained_intra_pred_flag = b.u1()
  pps.redundant_pic_cnt_present_flag = b.u1()

  if b.more_rbsp_data() {
    pps.transform_8x8_mode_flag = b.u1()
    pps.pic_scaling_matrix_present_flag = b.u1()
    if pps.pic_scaling_matrix_present_flag != 0 {
      for i in 0..6 + 2 * pps.transform_8x8_mode_flag {
        pps.pic_scaling_list_present_flag[i] = b.u1()
        if pps.pic_scaling_list_present_flag[i] != 0 {
          if i < 6 {
            // Make it a slice to get a dynamic array
            b.read_scaling_list(mut pps.scaling_list_4x4[i][0..], 16, mut &pps.use_default_scaling_matrix_4x4_flag[i])
          } else {
            b.read_scaling_list(mut pps.scaling_list_8x8[i - 6][0..], 64, mut &pps.use_default_scaling_matrix_8x8_flag[i - 6])
          }
        }
      }
    }
    pps.second_chroma_qp_index_offset = b.se()
  }
  b.read_rbsp_trailing_bits()
}

pub fn (sh SliceHeader) is_slice_type(cmp_type SH_SLICE_TYPE) bool {
  mut mslice_type := sh.slice_type
  mut mcmp_type := u32(cmp_type)
  if mslice_type >= 5 { mslice_type -= 5 }
  if mcmp_type >= 5 { mcmp_type -= 5 }
  if mslice_type == mcmp_type { return true }
  return false
}

pub fn (mut sh SliceHeader) read_ref_pic_list_reordering(mut b Bitstream) {
  if !sh.is_slice_type(SH_SLICE_TYPE.i) && !sh.is_slice_type(SH_SLICE_TYPE.si) {
    sh.rplr.ref_pic_list_reordering_flag_l0 = b.u1()
    if sh.rplr.ref_pic_list_reordering_flag_l0 != 0 {
      mut n := -1
      for {
        n++
        if n >= 64 {
          break
        }
        sh.rplr.reorder_l0.reordering_of_pic_nums_idc[n] = b.ue()
        if sh.rplr.reorder_l0.reordering_of_pic_nums_idc[n] == 0
        || sh.rplr.reorder_l0.reordering_of_pic_nums_idc[n] == 1 {
          sh.rplr.reorder_l0.abs_diff_pic_num_minus1[n] = b.ue()
        } else if sh.rplr.reorder_l0.reordering_of_pic_nums_idc[n] == 2 {
          sh.rplr.reorder_l0.long_term_pic_num[n] = b.ue()
        }
 
        if sh.rplr.reorder_l0.reordering_of_pic_nums_idc[n] == 3 || b.eof() {
          break
        }
      }
    }
  }

  if sh.is_slice_type(SH_SLICE_TYPE.b) {
    sh.rplr.ref_pic_list_reordering_flag_l1 = b.u1()
    if sh.rplr.ref_pic_list_reordering_flag_l1 != 0 {
      mut n := -1
      for {
        n++
        if n >= 64 {
          break
        }
        sh.rplr.reorder_l1.reordering_of_pic_nums_idc[n] = b.ue()
        if sh.rplr.reorder_l1.reordering_of_pic_nums_idc[n] == 0
        || sh.rplr.reorder_l1.reordering_of_pic_nums_idc[n] == 1 {
          sh.rplr.reorder_l1.abs_diff_pic_num_minus1[n] = b.ue()
        } else if sh.rplr.reorder_l1.reordering_of_pic_nums_idc[n] == 2 {
          sh.rplr.reorder_l1.long_term_pic_num[n] = b.ue()
        }
        if sh.rplr.reorder_l1.reordering_of_pic_nums_idc[n] == 3 || b.eof() {
          break
        }
      }
    }
  }
}

pub fn (mut sh SliceHeader) read_pred_weight_table(sps &SequenceParameterSet, pps &PictureParameterSet, mut b Bitstream) {
  sh.pwt.luma_log2_weight_denom = b.ue()
  if sps.chroma_format_idc != 0 {
    sh.pwt.chroma_log2_weight_denom = b.ue()
  }
  for i in 0..pps.num_ref_idx_l0_active_minus1 {    
    sh.pwt.luma_weight_l0_flag[i] = b.u1()
    if sh.pwt.luma_weight_l0_flag[i] != 0 {
      sh.pwt.luma_weight_l0[i] = b.se()
      sh.pwt.luma_offset_l0[i] = b.se()
    }
    if sps.chroma_format_idc != 0 {
      sh.pwt.chroma_weight_l0_flag[i] = b.u1()
      if sh.pwt.chroma_weight_l0_flag[i] != 0 {
        assert sh.pwt.chroma_weight_l0[i].len > 1
        for j in 0..2 {
          assert sh.pwt.chroma_weight_l0[i].len > 1
          sh.pwt.chroma_weight_l0[i][j] = b.se()
          sh.pwt.chroma_offset_l0[i][j] = b.se()
        }
      }
    }
  }

  if sh.is_slice_type(SH_SLICE_TYPE.b) {
    for i in 0..pps.num_ref_idx_l1_active_minus1 {
      sh.pwt.luma_weight_l1_flag[i] = b.u1()
      if sh.pwt.luma_weight_l1_flag[i] != 0 {
        sh.pwt.luma_weight_l1[i] = b.se()
        sh.pwt.luma_offset_l1[i] = b.se()
      }
      if sps.chroma_format_idc != 0 {
        sh.pwt.chroma_weight_l1_flag[i] = b.u1()
        if sh.pwt.chroma_weight_l1_flag[i] != 0 {
          assert sh.pwt.chroma_weight_l1[i].len > 1
          for j in 0..2 {
            assert sh.pwt.chroma_weight_l1[i].len > 1
            sh.pwt.chroma_weight_l1[i][j] = b.se()
            sh.pwt.chroma_offset_l1[i][j] = b.se()
          }
        }
      }
    }
  }
}

pub fn (mut sh SliceHeader) read_dec_ref_pic_marking(nal &NetworkAbstractionLayerHeader, mut b Bitstream) {
  if nal.type == NAL_UNIT_TYPE.coded_slice_idr {
    sh.drpm.no_output_of_prior_pics_flag = b.u1()
    sh.drpm.long_term_reference_flag = b.u1()
  } else {
    sh.drpm.adaptive_ref_pic_marking_mode_flag = b.u1()
    if sh.drpm.adaptive_ref_pic_marking_mode_flag != 0 {
      mut n := -1
      for {
        n++
        sh.drpm.memory_management_control_operation[n] = b.ue()
        if sh.drpm.memory_management_control_operation[n] == 1
        || sh.drpm.memory_management_control_operation[n] == 3 {
          sh.drpm.difference_of_pic_nums_minus1[n] = b.ue()
        }
        if sh.drpm.memory_management_control_operation[n] == 2 {
          sh.drpm.long_term_pic_num[n] = b.ue()
        }
        if sh.drpm.memory_management_control_operation[n] == 3
        || sh.drpm.memory_management_control_operation[n] == 6 {
          sh.drpm.long_term_frame_idx[n] = b.ue()
        }
        if sh.drpm.memory_management_control_operation[n] == 4 {
          sh.drpm.max_long_term_frame_idx_plus1[n] = b.ue()
        }
        
        if sh.drpm.memory_management_control_operation[n] == 0 || b.eof() {
          break
        }
      }
    }
  }
}

pub fn (mut sh SliceHeader) read_slice_header(nal &NetworkAbstractionLayerHeader, pps_array []PictureParameterSet, sps_array []SequenceParameterSet, mut b Bitstream) {
  sh.first_mb_in_slice = b.ue()
  sh.slice_type = b.ue()
  sh.pic_parameter_set_id = b.ue()

  pps := pps_array[sh.pic_parameter_set_id]
  sps := sps_array[pps.seq_parameter_set_id]
  
  sh.frame_num = b.u(u32(sps.log2_max_frame_num_minus4) + 4)

  if sps.frame_mbs_only_flag == 0 {
    sh.field_pic_flag = b.u1()
    if sh.field_pic_flag != 0 {
      sh.bottom_field_flag = b.u1()
    }
  }
  if nal.type == NAL_UNIT_TYPE.coded_slice_idr {
    sh.idr_pic_id = b.ue()
  }
  if sps.pic_order_cnt_type == 0 {
    sh.pic_order_cnt_lsb = b.u(u32(sps.log2_max_pic_order_cnt_lsb_minus4) + 4)
    if pps.pic_order_present_flag != 0 && sh.field_pic_flag == 0 {
      sh.delta_pic_order_cnt_bottom = b.se()
    }
  }
  if sps.pic_order_cnt_type == 1 && sps.delta_pic_order_always_zero_flag == 0 {
    sh.delta_pic_order_cnt[0] = b.se()
    if pps.pic_order_present_flag != 0 && sh.field_pic_flag == 0 {
      sh.delta_pic_order_cnt[1] = b.se()
    }
  }
  if pps.redundant_pic_cnt_present_flag != 0 {
    sh.redundant_pic_cnt = b.ue()
  }
  if sh.is_slice_type(SH_SLICE_TYPE.b) {
    sh.direct_spatial_mv_pred_flag = b.u1()
  }
  if sh.is_slice_type(SH_SLICE_TYPE.p)
  || sh.is_slice_type(SH_SLICE_TYPE.sp)
  || sh.is_slice_type(SH_SLICE_TYPE.b) {
    sh.num_ref_idx_active_override_flag = b.u1()
    if sh.num_ref_idx_active_override_flag != 0 {
      sh.num_ref_idx_l0_active_minus1 = b.ue()
      if sh.is_slice_type(SH_SLICE_TYPE.b) {
        sh.num_ref_idx_l1_active_minus1 = b.ue()
      }
    }
  }
  sh.read_ref_pic_list_reordering(mut b)
  if pps.weighted_pred_flag != 0
  && ((sh.is_slice_type(SH_SLICE_TYPE.p) || sh.is_slice_type(SH_SLICE_TYPE.sp))
  || (pps.weighted_bipred_idc == 1 && sh.is_slice_type(SH_SLICE_TYPE.b))) {
    sh.read_pred_weight_table(sps, pps, mut b)
  }
  if nal.idc != NAL_REF_IDC.priority_disposable {
    sh.read_dec_ref_pic_marking(nal, mut b)
  }
  if pps.entropy_coding_mode_flag != 0 && !sh.is_slice_type(SH_SLICE_TYPE.i) && !sh.is_slice_type(SH_SLICE_TYPE.si) {
    sh.cabac_init_idc = b.ue()
  }
  sh.slice_qp_delta = b.se()
  if sh.is_slice_type(SH_SLICE_TYPE.sp) || sh.is_slice_type(SH_SLICE_TYPE.si) {
    if sh.is_slice_type(SH_SLICE_TYPE.sp) {
      sh.sp_for_switch_flag = b.u1()
    }
    sh.slice_qs_delta = b.se()
  }
  if pps.deblocking_filter_control_present_flag != 0 {
    sh.disable_deblocking_filter_idc = b.ue()
    if sh.disable_deblocking_filter_idc != 1 {
      sh.slice_alpha_c0_offset_div2 = b.se()
      sh.slice_beta_offset_div2 = b.se()
    }
  }
  if pps.num_slice_groups_minus1 > 0
  && pps.slice_group_map_type >= 3
  && pps.slice_group_map_type <= 5 {
    v := intlog2(i32(pps.pic_size_in_map_units_minus1 + pps.slice_group_change_rate_minus1) + 1)
    sh.slice_group_change_cycle = b.u(u32(v))
  }
}


