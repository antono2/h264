# H.264 bitstream parser for V

This module parses the H.264/AVC syntax needed by Vulkan Video applications:
NAL headers, sequence and picture parameter sets (SPS/PPS), video usability
information (VUI), and slice headers.

It is **not** a video decoder: it does not perform entropy decoding, motion
compensation, inverse transforms, or produce pixels. The
[`v_vulkan_video`](https://github.com/antono2/v_vulkan_video) player uses these
parsed structures to prepare hardware decode operations.

## Install

```bash
v install https://github.com/antono2/h264
```

## Minimal example

```v
import h264

mut stream := h264.Bitstream{}
stream.init([u8(0x67)]) // forbidden_zero_bit=0, nal_ref_idc=3, type=SPS

mut header := h264.NetworkAbstractionLayerHeader{}
header.read_nal_header(mut stream)
assert header.type == .sequence_parameter_set
```

`Bitstream.init()` expects raw RBSP/NAL bytes in memory. Container extraction,
length prefixes or Annex-B start codes, emulation-prevention removal, frame
reordering, and decoded-picture management remain the caller's responsibility.

## Safety and supported scope

The parser exposes low-level syntax structures rather than a defensive media
API. Some invalid syntax is rejected with assertions, while truncated fields
can read as zero at end of input. Validate untrusted container lengths and NAL
boundaries before parsing them.

The implementation covers the syntax exercised by the Vulkan Video H.264
player and is not yet a claim of complete support for every profile, extension,
bit depth, chroma format, or interlaced stream.

## Tests

The included tests are software-only and require no GPU:

```bash
v test .
```

They cover RBSP look-ahead behavior and representative High-profile SPS/PPS
parsing used by the player.
