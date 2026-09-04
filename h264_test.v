module h264

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
