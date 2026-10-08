// Verifies HRD entry counts and the cursor position of following syntax fields.
module h264

import encoding.hex

fn test_hrd_reads_all_cpb_entries() {
	// H.264 E.1.2 encodes count minus one. Fixtures contain 1, 2 and 32 CPBs,
	// distinct entry values, four delay widths and an 0xa5 cursor sentinel.
	for count, encoded in {
		1:  '9a30208a63a294'
		2:  '468c08220102531d14a0'
		32: '04068c08220102507e303e9c1e8403c89076140ea2c1c860388d06e1c0da3c1a8200d2110660903284c188280c21505e0b02e85c168300b2190560d02a86c148380a21d04e0f02687c1281002482104604408a531d14a0'
	} {
		mut stream := Bitstream{}
		stream.init(hex.decode(encoded)!)
		mut sps := SequenceParameterSet{}
		sps.read_hrd_parameters(mut stream)
		assert sps.hrd.cpb_cnt_minus1 == u32(count - 1)
		assert sps.hrd.bit_rate_scale == 3
		assert sps.hrd.cpb_size_scale == 4
		for i in 0 .. count {
			assert sps.hrd.bit_rate_value_minus1[i] == u32(i + 2)
			assert sps.hrd.cpb_size_value_minus1[i] == u32(64 - i)
			assert sps.hrd.cbr_flag[i] == u32(i % 2)
		}
		assert sps.hrd.initial_cpb_removal_delay_length_minus1 == 5
		assert sps.hrd.cpb_removal_delay_length_minus1 == 6
		assert sps.hrd.dpb_output_delay_length_minus1 == 7
		assert sps.hrd.time_offset_length == 8
		assert stream.u(8) == 0xa5
	}
}
