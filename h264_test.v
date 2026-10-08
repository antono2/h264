// Exercises bit reads, NAL headers and parameter-set parsing without a decoder or GPU.
// Includes scaling-list regressions and checks for EOF and cursor-preserving look-ahead.
module h264

import encoding.hex

fn test_custom_sps_scaling_lists_are_stored() {
	// FRExt_MMCO4_Sony_B conformance SPS, excluding its NAL header.
	data :=
		hex.decode('64001fad9464763b8ac4444a323b1dc5622225191d8ee2b11114222b373669a844566e6cd35088acdcd9a69444cd1b9bc57c9f93f9bf27c9e4e4cd251a4689c9ebe4fd7f27ebe4f5c9a906c694160964')!
	mut stream := Bitstream{}
	stream.init(data)
	mut sps := SequenceParameterSet{}
	sps.read_sps(mut stream)
	assert sps.seq_parameter_set_id == 0
	assert sps.scaling_list_4x4[0][..4] == [i32(6), 12, 12, 19]
	assert sps.scaling_list_8x8[0][..4] == [i32(6), 10, 10, 13]
	for i in 0 .. 6 {
		assert sps.scaling_list_4x4[i][15] != 0
	}
	for i in 0 .. 2 {
		assert sps.scaling_list_8x8[i][63] != 0
	}
	// Check the complete persisted record, not sizeof(variable), which current
	// V3 can evaluate as a heap-promoted pointer's size.
	mut stored := []u8{}
	stored << unsafe { byteptr(&sps).vbytes(int(sizeof(SequenceParameterSet))) }
	assert stored.len == int(sizeof(SequenceParameterSet))
	record := unsafe { &SequenceParameterSet(stored.data) }
	assert record.scaling_list_4x4[0][..4] == [i32(6), 12, 12, 19]
	assert record.scaling_list_8x8[0][..4] == [i32(6), 10, 10, 13]
}

fn test_custom_pps_scaling_list_is_stored() {
	mut stream := Bitstream{}
	stream.init(hex.decode('ce3c7fffe0c0')!)
	mut pps := PictureParameterSet{}
	pps.read_pps(mut stream)
	assert pps.pic_scaling_matrix_present_flag == 1
	for value in pps.scaling_list_4x4[0] {
		assert value == 8
	}
}

fn test_scaling_list_default_flag_points_to_value_not_pointer() {
	mut stream := Bitstream{}
	// signed Exp-Golomb -8 makes the first nextScale zero.
	stream.init([u8(0x08), 0x80])
	mut list := []i32{len: 16}
	mut use_default := u32(0)
	stream.read_scaling_list(mut list, 16, &use_default)
	assert use_default == 1
	for value in list {
		assert value == 8
	}
}

fn test_fixed_width_reads_cross_byte_boundaries() {
	mut stream := Bitstream{}
	stream.init([u8(0xa5), 0xc0])
	assert stream.u(8) == 0xa5
	assert stream.byte_aligned()
	assert stream.u(2) == 0x3
	assert !stream.byte_aligned()
}

fn test_truncated_reads_are_zero_filled() {
	mut stream := Bitstream{}
	stream.init([]u8{})
	assert stream.u(16) == 0
	assert stream.eof()
}

fn test_unsigned_exp_golomb_values() {
	for sample in [
		[u8(0b10000000)],
		[u8(0b01000000)],
		[u8(0b01100000)],
		[u8(0b00100000)],
	] {
		mut stream := Bitstream{}
		stream.init(sample)
		expected := u32(match sample[0] {
			0b10000000 { 0 }
			0b01000000 { 1 }
			0b01100000 { 2 }
			else { 3 }
		})
		assert stream.ue() == expected
	}
}

fn test_signed_exp_golomb_values() {
	for sample, expected in {
		u8(0b10000000): 0
		u8(0b01000000): 1
		u8(0b01100000): -1
		u8(0b00100000): 2
		u8(0b00101000): -2
	} {
		mut stream := Bitstream{}
		stream.init([sample])
		assert stream.se() == expected
	}
}

fn test_read_nal_header() {
	mut stream := Bitstream{}
	stream.init([u8(0x67)])
	mut header := NetworkAbstractionLayerHeader{}
	header.read_nal_header(mut stream)
	assert header.idc == .priority_highest
	assert header.type == .sps
}

fn test_intlog2_uses_minimum_required_bits() {
	assert intlog2(0) == 0
	assert intlog2(1) == 0
	assert intlog2(2) == 1
	assert intlog2(3) == 2
	assert intlog2(4) == 2
	assert intlog2(5) == 3
}

fn test_more_rbsp_data_preserves_bit_position() {
	mut stream := Bitstream{}
	stream.init([u8(0b10100000)])
	stream.u(2)
	start := stream.p
	bits_left := stream.bits_left

	assert !stream.more_rbsp_data()
	assert stream.p == start
	assert stream.bits_left == bits_left
}

fn test_parse_sample_high_profile_pps() {
	mut stream := Bitstream{}
	// The MP4 sample's PPS RBSP, after its 0x68 NAL header.
	stream.init([u8(0xee), 0x0d, 0x8b])
	mut pps := PictureParameterSet{}
	pps.read_pps(mut stream)
	assert pps.pic_parameter_set_id == 0
	assert pps.seq_parameter_set_id == 0
	mut stored_pps := []u8{}
	stored_pps.ensure_cap(int(sizeof(PictureParameterSet)))
	stored_pps << unsafe { byteptr(&pps).vbytes(int(sizeof(PictureParameterSet))) }
	assert stored_pps.len == int(sizeof(PictureParameterSet))
}

fn test_persist_sample_parameter_sets() {
	mut sps_stream := Bitstream{}
	sps_stream.init([u8(0x64), 0x00, 0x28, 0xac, 0xb4, 0x03, 0xc0, 0x11, 0x3f, 0x2c, 0xd4, 0x04,
		0x04, 0x04, 0x1e, 0x2c, 0x5d, 0x40])
	mut sps := SequenceParameterSet{}
	sps.read_sps(mut sps_stream)
	assert sps.profile_idc == 100
	assert sps.level_idc == 40
	assert sps.seq_parameter_set_id == 0
	sps.seq_parameter_set_id = 7
	mut stored_sps := []u8{}
	stored_sps << unsafe { byteptr(&sps).vbytes(int(sizeof(SequenceParameterSet))) }

	mut pps_stream := Bitstream{}
	pps_stream.init([u8(0xee), 0x0d, 0x8b])
	mut pps := PictureParameterSet{}
	pps.read_pps(mut pps_stream)
	pps.pic_parameter_set_id = 11
	pps.seq_parameter_set_id = 7
	mut stored_pps := []u8{}
	stored_pps << unsafe { byteptr(&pps).vbytes(int(sizeof(PictureParameterSet))) }

	assert stored_sps.len == int(sizeof(SequenceParameterSet))
	assert stored_pps.len == int(sizeof(PictureParameterSet))
	stored_sps_record := unsafe { &SequenceParameterSet(stored_sps.data) }
	stored_pps_record := unsafe { &PictureParameterSet(stored_pps.data) }
	assert stored_sps_record.seq_parameter_set_id == 7
	assert stored_pps_record.pic_parameter_set_id == 11
	assert stored_pps_record.seq_parameter_set_id == 7
}
