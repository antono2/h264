# H.264 bitstream parser for V

[Project portfolio](https://oreskin.de/projects_en.php)

This module parses the H.264/AVC syntax needed by Vulkan Video applications:
NAL headers, sequence and picture parameter sets (SPS/PPS), video usability
information (VUI) and slice headers.

It is **not** a video decoder: it does not perform entropy decoding, motion compensation or inverse transforms, and it does not produce pixels. The
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
assert header.type == .sps
```

`Bitstream.init()` expects raw RBSP/NAL bytes in memory. Container extraction,
length prefixes or Annex-B start codes, emulation-prevention removal, frame
reordering and decoded-picture management remain the caller's responsibility.

## Parsing order and caller responsibilities

1. Extract one complete NAL unit from its container or Annex-B stream. Remove
   the container length or Annex-B prefix and remove emulation-prevention bytes
   from the payload before parsing it.
2. Initialize a `Bitstream` with the NAL header followed by the unescaped payload,
   then call `read_nal_header`. The cursor now points at the payload. If you
   already separated the header, initialize a payload-only stream for the
   corresponding parameter-set or slice parser instead.
3. For `.sps` and `.pps`, parse into a fresh `SequenceParameterSet` or
   `PictureParameterSet`. Preserve each result by its syntax ID; optional fields
   are not all reset when reusing an existing value.
4. For a supported coded slice, call `read_slice_header` with its NAL header and
   SPS/PPS arrays indexed by `seq_parameter_set_id` and `pic_parameter_set_id`.
   Every referenced entry must contain a parsed parameter set. IDs outside the
   supplied array lengths cause a panic; array bounds alone do not establish
   that an entry is populated.

`Bitstream.init` retains the supplied bytes without cloning them. Keep those
bytes available and unchanged while parsing, and give concurrent parses
separate cursors and output structures. Parameter-set readers consume trailing
RBSP bits; the slice-header reader leaves the remaining coded slice data for
the decoder. Returned fields retain their H.264 syntax encodings rather than
being normalized into dimensions, frame rates or decoder operations.

## Safety and supported scope

The parser exposes low-level syntax structures rather than a defensive media
API. Some invalid syntax is rejected with assertions, while truncated fields
can read as zero at end of input. Validate untrusted container lengths and NAL
boundaries before parsing them.

The implementation covers the syntax exercised by the Vulkan Video H.264
player and is not yet a claim of complete support for every profile, extension,
bit depth, chroma format or interlaced stream.

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
NAL headers, RBSP look-ahead, integer bit widths and representative
High-profile SPS/PPS parsing used by the player. HRD buffering syntax tests cover
one, two and the maximum 32 CPB entries, including alignment of the following
delay fields. Counts beyond 32 entries are rejected by an assertion.

CI also runs strict V3 with TinyCC on Linux, using a pinned compiler snapshot.
The regression samples cover custom SPS/PPS scaling lists and the default-list
flag. Scaling lists are copied explicitly into their fixed storage arrays.

Callers of `Bitstream.read_scaling_list` should pass the flag as `&flag`, not
`mut &flag`. The pointer itself is not reassigned; only its pointed-to value is
written. This avoids a pointer-to-pointer mismatch in V3-generated C.
