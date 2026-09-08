module h264

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
	stored_pps.ensure_cap(int(sizeof(pps)))
	stored_pps << unsafe { byteptr(&pps).vbytes(int(sizeof(pps))) }
	assert stored_pps.len == int(sizeof(pps))
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
	mut stored_sps := []u8{}
	stored_sps << unsafe { byteptr(&sps).vbytes(int(sizeof(sps))) }

	mut pps_stream := Bitstream{}
	pps_stream.init([u8(0xee), 0x0d, 0x8b])
	mut pps := PictureParameterSet{}
	pps.read_pps(mut pps_stream)
	mut stored_pps := []u8{}
	stored_pps << unsafe { byteptr(&pps).vbytes(int(sizeof(pps))) }

	assert stored_sps.len == int(sizeof(sps))
	assert stored_pps.len == int(sizeof(pps))
}
