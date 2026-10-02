# H.264 bitstream parser for V

[Project portfolio](https://oreskin.de/projects_en.php)

This module parses the H.264/AVC syntax needed by Vulkan Video applications:
NAL headers, sequence and picture parameter sets (SPS/PPS), video usability
information (VUI), and slice headers.

It is **not** a video decoder: it does not perform entropy decoding, motion
compensation, inverse transforms, or produce pixels. The
[`v_vulkan_video`](https://github.com/antono2/v_vulkan_video) player uses these
parsed structures to prepare hardware decode operations.

## Install

```bash
v install antono2.h264
```

## Platform support

The parser is implemented entirely in V and has no native library or GPU
dependency. CI runs the complete test suite on Linux and on Windows Server 2022
with MSVC. Windows 10/11 x64 users can install and import the module with the
same `v install` command shown above.

## Minimal example

```v
import antono2.h264

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

The included tests are software-only and require no GPU. Run them with the
default compiler, or select MSVC on Windows:

```bash
v test .
```

```powershell
v -cc msvc test .
```

They cover fixed-width and Exp-Golomb bit reading, truncated input behavior,
NAL headers, RBSP look-ahead, integer bit widths, and representative
High-profile SPS/PPS parsing used by the player.

CI also runs strict V3 with TinyCC on Linux, using a pinned compiler snapshot.
The regression samples cover custom SPS/PPS scaling lists and the default-list
flag. Scaling lists are copied explicitly into their fixed storage arrays.

Callers of `Bitstream.read_scaling_list` should pass the flag as `&flag`, not
`mut &flag`. The pointer itself is not reassigned; only its pointed-to value is
written. This avoids a pointer-to-pointer mismatch in V3-generated C.
