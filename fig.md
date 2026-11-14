# `.fig`

```toml
[fig]
media_type = "text/vnd.fig+markdown"
field_type = "application/vnd.fig+raw"
extension  = ".fig"
version    = "1.0.0"

[pipeline]
input   = "biomedical telemetry --multi-dimensional"
process = "z-order { finite-field GF(2^8)^w, scale-invariant, gif-time }"
output  = "graph.coord --fileformat .fig { SUPERSEDES .dicom, .cif }"
```

## 1. Axioms

1. **Superset.** Every CommonMark `.md` or OpenGigantic `.og` file is a valid `.fig` link graph. Entry is `::index`. `:requires` loads dependencies. `::warning` is never skipped. `scope` gates visibility. Graphs compose across HTTP and IPFS `uri#id` edges.
2. **Headerless field.** Terminal sample streams carry zero header bytes, magic signatures, dimension tags, or palettes. Byte `0` is the least-significant byte of the first record `le(code, c) || word`, ordered by strictly increasing integer Morton `code` (`z_0 < z_1 < ... < z_{n-1}`).
3. **Content address & chain identity.** Field and figure identities are CIDv1 raw SHA-256 (`0x01 || 0x55 || 0x12 || 0x20 || sha256(bytes)`, multibase `base32lower`). Every HTTP `GET` body must match the `cid` declared on its node; the 256-bit digest maps bijectively to `uint256` `tokenId` alongside the `.fig` `tokenURI`.
4. **Z-order isomorphism.** Over `d >= 1` spatial axes at bit depth `b >= 1` (`d * b <= 256`, `c = ceil(d * b / 8)` bytes), Morton interleaving `Z : (GF(2)^b)^d -> GF(2)^(d*b)` places axis `0` at bit `0` of each `d`-bit stratum: `z(x) = sum_{k=0..b-1} sum_{i=0..d-1} bit_k(x_i) * 2^(k*d + i)`.
5. **Scale invariance.** Dropping `k` bit-planes shifts `z |-> z >> (k * d)` and scales the rational metre quantum `q |-> q * 2^k`. Every `2^k`-aligned block spans `2^(k*d)` contiguous codes. At `q = 10^-12 m` (`1 pm`) and `b = 48`, one unsigned lattice addresses covalent bonds (`154 pm`) through whole-body anatomy (`1.8 m`) and coarsens at `k = 30` to millimetre voxel grids (`2^30 pm ≈ 1.074 mm`).
6. **Finite-field raw words.** Every `w`-byte sample (`w >= 1`) is a raw binary vector in `(GF(2^8))^w` modulo `x^8 + x^4 + x^3 + x + 1` (`0x11b`): vector addition is bitwise XOR and scalar multiplication reduces modulo `0x11b`. Physical calibration applies a rational affine reading `slope * leInt(word) + intercept` to the same raw bytes without altering field operations.
7. **GIF temporal primitives.** Frames borrow GIF89a graphic control and Netscape 2.0 looping: `delay` in rational `time-num / time-den` seconds (default `1/100 s`), `disposal` (`0` normalized to `1` keep, `2` restore graphic domain to `background`, `3` restore overwritten cells to prior canvas), optional `transparent` `w`-byte word, `user-input` (`1` halts autonomous progression), and `loop` (`N > 0` plays `N` times; `0` unbounded under an explicit caller cap). Record absence, `background`, and `transparent` stay distinct.
8. **Supersession (`.dicom`, `.cif`).** Typed `embedding` nodes project `.dicom` patient space (`IPP`, `IOP`, spacing, rescale) and `.cif` crystal space (`fract mod 1`, unit-cell metric) onto the lattice by round-half-away quantization `roundHalfAway((x_m - origin) / q)`, rejecting out-of-range coordinates, quantization collisions, and patient-identifying headers.

## 2. Ontology & Coordinate Algebra (`lean`)

| Node | Role | Edge | Role |
|---|---|---|---|
| `figure` | Root (`loop`, `background`, `time-num`, `time-den`) | `:contains` | Child scope (`lattice`, `alphabet`, `frame`, `field`, `embedding`) |
| `lattice` | Finite grid (`dims`, `bits`, `axes`, `quantum-num`, `quantum-den`, `origin`) | `:requires` | Premise, lattice, alphabet, or guard dependency |
| `alphabet` | Finite field `(GF(2^8))^w` (`poly=0x11b`, `value-bytes`) | `:precedes` | Strict frame order (`seq = 0 .. n-1`) |
| `frame` | GIF control (`seq`, `delay`, `disposal`, `transparent`, `user-input`, `bounds`) | `:links` | HTTP or IPFS locator verified against `cid` |
| `field` | Headerless stream (`cid`, inline `hex` or remote raw bytes) | `:coarsens` | Renormalized lattice view (`levels=k`, `z >> k*d`, `q * 2^k`) |
| `embedding` | Coordinate projection (`kind=dicom\|cif`) | `:supersedes` | Replaced `.dicom` or `.cif` coordinate payload |
| `invariant` | Safety, canonicality, or content-address gate | | |
| `meta` | `__view`, `__changelog` | | |

```lean
inductive NodeType | figure | lattice | alphabet | frame | field | embedding | invariant | «meta»
  deriving DecidableEq, Repr
inductive EdgeType | contains | requires | precedes | links | coarsens | supersedes
  deriving DecidableEq, Repr
inductive Disposal | unspecified | keep | background | previous
  deriving DecidableEq, Repr
inductive EmbedKind | dicom | cif
  deriving DecidableEq, Repr

structure Q where
  num : Int
  den : Int
  deriving DecidableEq, Repr

namespace Q
def norm (q : Q) : Q :=
  let g : Int := Int.gcd q.num q.den
  let s : Int := if q.den < 0 then -1 else 1
  ⟨s * (q.num / g), s * (q.den / g)⟩
def add (a b : Q) : Q := norm ⟨a.num * b.den + b.num * a.den, a.den * b.den⟩
def sub (a b : Q) : Q := norm ⟨a.num * b.den - b.num * a.den, a.den * b.den⟩
def mul (a b : Q) : Q := norm ⟨a.num * b.num, a.den * b.den⟩
def div (a b : Q) : Q := norm ⟨a.num * b.den, a.den * b.num⟩
def mod1 (q : Q) : Q :=
  let n := norm q
  let r := n.num % n.den
  let r' := if r < 0 then r + n.den else r
  norm ⟨r', n.den⟩
def roundHalfAway (q : Q) : Int :=
  let n := norm q
  let s : Int := if n.num < 0 then -1 else 1
  let a : Int := Int.natAbs n.num
  let quot := a / n.den
  let rem  := a % n.den
  s * (if 2 * rem ≥ n.den then quot + 1 else quot)
def quantize (val origin quantum : Q) : Int :=
  roundHalfAway (div (sub val origin) quantum)
end Q

def codeBytes (dims bits : Nat) : Nat := (dims * bits + 7) / 8

def mortonEncode : Nat → List Nat → Nat
  | 0, _ => 0
  | b + 1, xs =>
      let plane := xs.foldr (fun x acc => x % 2 + 2 * acc) 0
      plane + (2 ^ xs.length) * mortonEncode b (xs.map (· / 2))

def extractPlane : Nat → Nat → List Nat × Nat
  | 0, z => ([], z)
  | d + 1, z =>
      let (rest, z') := extractPlane d (z / 2)
      ((z % 2) :: rest, z')

def mortonDecode : Nat → Nat → Nat → List Nat
  | _, d, 0 => List.replicate d 0
  | z, d, b + 1 =>
      let (lo, z') := extractPlane d z
      List.zipWith (fun bit hi => bit + 2 * hi) lo (mortonDecode z' d b)

def coarsen : Nat → Nat → Nat → Nat
  | z, _, 0 => z
  | z, d, k + 1 => coarsen (z / (2 ^ d)) d k

def gfAdd (a b : Nat) : Nat := (a % 256) ^^^ (b % 256)

def gfMul (a b : Nat) : Nat :=
  let rec go : Nat → Nat → Nat → Nat → Nat
    | 0, _, _, p => p
    | i + 1, a, b, p =>
        let p := if b % 2 == 1 then p ^^^ a else p
        let hi := a / 128
        let a := (a * 2) % 256
        let a := if hi == 1 then a ^^^ 0x1b else a
        go i a (b / 2) p
  go 8 (a % 256) (b % 256) 0

def gfVecAdd (xs ys : List Nat) : List Nat := List.zipWith gfAdd xs ys
def gfScale (c : Nat) (xs : List Nat) : List Nat := xs.map (gfMul c)

def leBytes : Nat → Nat → List Nat
  | _, 0 => []
  | n, w + 1 => (n % 256) :: leBytes (n / 256) w

def leInt : List Nat → Int
  | [] => 0
  | b :: bs => Int.ofNat (b % 256) + 256 * leInt bs

def fitsBytes : Nat → Nat → Bool
  | n, 0 => n == 0
  | n, w + 1 => fitsBytes (n / 256) w

def rescale (word : List Nat) (slope intercept : Q) : Q :=
  Q.add (Q.mul slope ⟨leInt word, 1⟩) intercept

structure Record where
  code : Nat
  word : List Nat
  deriving DecidableEq, Repr

def strictlyIncreasing : List Nat → Bool
  | [] => true
  | [_] => true
  | a :: b :: rest => a < b && strictlyIncreasing (b :: rest)

def serializeField (codeW wordW : Nat) (rs : List Record) : Option (List Nat) :=
  if strictlyIncreasing (rs.map Record.code) &&
     rs.all (fun r => fitsBytes r.code codeW && r.word.length == wordW && wordW > 0) then
    some (rs.foldl (fun acc r => acc ++ leBytes r.code codeW ++ r.word) [])
  else none

structure Frame where
  samples     : List Record
  delay       : Nat
  disposal    : Disposal
  transparent : Option (List Nat)
  userInput   : Bool
  bounds      : Option (List Nat)
  deriving DecidableEq, Repr

structure Moment where
  visible : List Record
  rested  : List Record
  deriving DecidableEq, Repr

def canvasGet (c : List Record) (bg : List Nat) (z : Nat) : List Nat :=
  match c.find? (fun r => r.code == z) with
  | some r => r.word
  | none   => bg

def canvasPut (c : List Record) (bg : List Nat) (z : Nat) (w : List Nat) : List Record :=
  let rest := c.filter (fun r => r.code != z)
  if w == bg then rest else rest ++ [⟨z, w⟩]

def stepFrame (canvas : List Record) (bg : List Nat) (f : Frame) : Option Moment :=
  if f.userInput then none
  else
    let drawn := f.samples.filter (fun r => f.transparent != some r.word)
    let vis := drawn.foldl (fun c r => canvasPut c bg r.code r.word) canvas
    let domain := match f.bounds with
      | some zs => zs
      | none    => f.samples.map Record.code
    let rest := match f.disposal with
      | .unspecified | .keep => vis
      | .background => domain.foldl (fun c z => canvasPut c bg z bg) vis
      | .previous   => drawn.foldl (fun c r => canvasPut c bg r.code (canvasGet canvas bg r.code)) vis
    some ⟨vis, rest⟩

def playFrames (bg : List Nat) : List Record → List Frame → Option (List Moment)
  | _, [] => some []
  | c, f :: fs =>
      match stepFrame c bg f with
      | none => none
      | some m =>
          match playFrames bg m.rested fs with
          | none => none
          | some ms => some (m :: ms)

def fitsAxis : Int → Nat → Bool
  | x, 0 => x == 0
  | x, b + 1 => 0 ≤ x && fitsAxis (x / 2) b

def dicomVoxel (col row slc : Int) (ipp cDir rDir sDir : List Q)
    (cSp rSp sSp : Q) (origin : List Q) (quantum : Q) (bits : Nat) : Option (List Nat) :=
  let coords := (List.range 3).map (fun a =>
    let mm := Q.add (ipp.getD a ⟨0, 1⟩)
      (Q.add (Q.mul ⟨col, 1⟩ (Q.mul cSp (cDir.getD a ⟨0, 1⟩)))
      (Q.add (Q.mul ⟨row, 1⟩ (Q.mul rSp (rDir.getD a ⟨0, 1⟩)))
             (Q.mul ⟨slc, 1⟩ (Q.mul sSp (sDir.getD a ⟨0, 1⟩)))))
    let m := Q.div mm ⟨1000, 1⟩
    Q.quantize m (origin.getD a ⟨0, 1⟩) quantum)
  if coords.all (fun x => fitsAxis x bits) then
    some (coords.map Int.toNat)
  else none

def cifOrtho (frac cellA : List Q) (origin : List Q) (quantum : Q) (bits : Nat) : Option (List Nat) :=
  let ang : Q := ⟨1, 10000000000⟩
  let coords := (List.range 3).map (fun a =>
    let m := Q.mul (Q.mul (Q.mod1 (frac.getD a ⟨0, 1⟩)) (cellA.getD a ⟨0, 1⟩)) ang
    Q.quantize m (origin.getD a ⟨0, 1⟩) quantum)
  if coords.all (fun x => fitsAxis x bits) then
    some (coords.map Int.toNat)
  else none

def injectiveCoords : List (List Nat) → Bool
  | [] => true
  | x :: xs => !(xs.contains x) && injectiveCoords xs

def roundtripAll (dims bits : Nat) : Bool :=
  (List.range (2 ^ (dims * bits))).all (fun z =>
    let cs := mortonDecode z dims bits
    mortonEncode bits cs == z && cs.all (fun x => fitsAxis (Int.ofNat x) bits))

theorem zorder_laws :
    roundtripAll 2 2 = true ∧
    roundtripAll 3 2 = true ∧
    mortonEncode 1 [1, 0, 0] = 1 ∧
    mortonEncode 1 [0, 1, 0] = 2 ∧
    mortonEncode 1 [0, 0, 1] = 4 ∧
    mortonEncode 2 [3, 2] = 13 ∧
    ([((2 : Nat), (2 : Nat)), (3, 2), (2, 3), (3, 3)].map (fun (x, y) => mortonEncode 2 [x, y])) =
      [12, 13, 14, 15] ∧
    coarsen (mortonEncode 4 [3, 2]) 2 1 = mortonEncode 3 [1, 1] ∧
    codeBytes 3 48 = 18 ∧
    mortonDecode (mortonEncode 48 [154, 0, 0]) 3 48 = [154, 0, 0] ∧
    mortonDecode (mortonEncode 48 [1800000000000, 0, 0]) 3 48 = [1800000000000, 0, 0] ∧
    coarsen (mortonEncode 48 [1800000000000, 0, 0]) 3 30 = mortonEncode 18 [1676, 0, 0] := by decide

theorem gf256_laws :
    gfMul 0x57 0x83 = 0xc1 ∧
    gfMul 0x57 0xbf = 0x01 ∧
    gfAdd 0x57 0x83 = 0xd4 ∧
    (0x57 + 0x83) % 256 = 0xda ∧
    gfMul (gfAdd 0x57 0x83) 0x13 = gfAdd (gfMul 0x57 0x13) (gfMul 0x83 0x13) ∧
    gfVecAdd [0x57, 0x01] [0x83, 0x01] = [0xd4, 0x00] ∧
    gfScale 0x57 [0x01, 0x83] = [0x57, 0xc1] ∧
    rescale [0x90, 0x01] ⟨1, 1⟩ ⟨-1024, 1⟩ = ⟨-624, 1⟩ := by decide

def f0 : Frame := ⟨[⟨0, [0x57]⟩, ⟨1, [0x83]⟩, ⟨2, [0xd4]⟩, ⟨3, [0xc1]⟩], 10, .keep, none, false, none⟩
def f1 : Frame := ⟨[⟨3, [0x00]⟩, ⟨12, [0x11]⟩], 5, .background, some [0x00], false, none⟩
def f1Prev : Frame := { f1 with disposal := .previous }
def fRegion : Frame := ⟨[⟨0, [0x11]⟩], 1, .background, none, false, some [0, 1]⟩

theorem stream_and_gif_laws :
    serializeField 1 1 f0.samples = some [0x00, 0x57, 0x01, 0x83, 0x02, 0xd4, 0x03, 0xc1] ∧
    serializeField 1 1 f1.samples = some [0x03, 0x00, 0x0c, 0x11] ∧
    serializeField 2 1 [⟨1, [0xaa]⟩, ⟨256, [0xbb]⟩] = some [0x01, 0x00, 0xaa, 0x00, 0x01, 0xbb] ∧
    serializeField 2 1 [⟨256, [0xbb]⟩, ⟨1, [0xaa]⟩] = none ∧
    serializeField 1 1 [⟨3, [0x01]⟩, ⟨3, [0x02]⟩] = none ∧
    stepFrame [] [0x00] { f0 with disposal := .unspecified } = stepFrame [] [0x00] f0 ∧
    (match playFrames [0x00] [] [f0, f1] with
     | some [_, m1] =>
         canvasGet m1.visible [0x00] 3 == [0xc1] &&
         canvasGet m1.visible [0x00] 12 == [0x11] &&
         canvasGet m1.rested [0x00] 3 == [0x00] &&
         canvasGet m1.rested [0x00] 12 == [0x00] &&
         canvasGet m1.rested [0x00] 2 == [0xd4]
     | _ => false) = true ∧
    (match playFrames [0x00] [] [f0, f1Prev] with
     | some [_, m1] =>
         canvasGet m1.rested [0x00] 3 == [0xc1] &&
         canvasGet m1.rested [0x00] 12 == [0x00]
     | _ => false) = true ∧
    (match playFrames [0x00] [] [f0, fRegion] with
     | some [_, m1] =>
         canvasGet m1.rested [0x00] 0 == [0x00] &&
         canvasGet m1.rested [0x00] 1 == [0x00] &&
         canvasGet m1.rested [0x00] 2 == [0xd4]
     | _ => false) = true ∧
    stepFrame [] [0x00] { f0 with userInput := true } = none := by decide

def z3 : List Q := [⟨0, 1⟩, ⟨0, 1⟩, ⟨0, 1⟩]
def ex : List Q := [⟨1, 1⟩, ⟨0, 1⟩, ⟨0, 1⟩]
def ey : List Q := [⟨0, 1⟩, ⟨1, 1⟩, ⟨0, 1⟩]
def ez : List Q := [⟨0, 1⟩, ⟨0, 1⟩, ⟨1, 1⟩]

theorem embed_laws :
    dicomVoxel 1 1 3 [⟨3, 2⟩, ⟨0, 1⟩, ⟨0, 1⟩] ex ey ez ⟨2, 1⟩ ⟨5, 1⟩ ⟨4, 1⟩ z3 ⟨1, 1000⟩ 16 =
      some [4, 5, 12] ∧
    dicomVoxel 0 0 0 [⟨-10, 1⟩, ⟨0, 1⟩, ⟨0, 1⟩] ex ey ez ⟨1, 1⟩ ⟨1, 1⟩ ⟨1, 1⟩ z3 ⟨1, 1000⟩ 16 =
      none ∧
    dicomVoxel 0 0 0 [⟨-10, 1⟩, ⟨0, 1⟩, ⟨0, 1⟩] ex ey ez ⟨1, 1⟩ ⟨1, 1⟩ ⟨1, 1⟩
      [⟨-1, 100⟩, ⟨0, 1⟩, ⟨0, 1⟩] ⟨1, 1000⟩ 16 = some [0, 0, 0] ∧
    injectiveCoords [[0, 0, 0], [0, 0, 0]] = false ∧
    cifOrtho [⟨1, 2⟩, ⟨1, 5⟩, ⟨-1, 5⟩] [⟨10, 1⟩, ⟨10, 1⟩, ⟨10, 1⟩] z3 ⟨1, 10000000000⟩ 8 =
      some [5, 2, 8] ∧
    cifOrtho [⟨1, 1⟩, ⟨1, 2⟩, ⟨0, 1⟩] [⟨5, 1⟩, ⟨5, 1⟩, ⟨5, 1⟩] z3 ⟨1, 10000000000⟩ 8 =
      some [0, 3, 0] ∧
    decide (87 * 87 - 7500 < 7500 - 86 * 86) = true := by decide
```

Non-orthogonal `.cif` cells (`(alpha, beta, gamma) /= (90, 90, 90)`) wrap fractional coordinates into `[0, 1)` and apply the upper-triangular PDB orthogonalization `(X, Y, Z)_A = (a*xf + b*cos(gamma)*yf + c*cos(beta)*zf, b*sin(gamma)*yf + c*(cos(alpha) - cos(beta)*cos(gamma))/sin(gamma)*zf, c*sqrt(1 - cos^2(alpha) - cos^2(beta) - cos^2(gamma) + 2*cos(alpha)*cos(beta)*cos(gamma))/sin(gamma)*zf)` in IEEE-754 binary64 before `roundHalfAway(((X, Y, Z) * 10^-10 - origin) / quantum)`; degenerate cells (`sin(gamma) = 0` or non-positive volume) are rejected. Chunked HTTP transport may attach an optional RFC 6962 SHA-256 Merkle root (`0x00 || chunk` leaves, `0x01 || left || right` parents); the canonical identity remains the whole-field CIDv1 raw digest.

## 3. Grammar (`~EBNF`)

```ebnf
fig          = markdown_raw | manifest , { node_block | markdown_prose } ;
manifest     = "::meta"   , NL , { kv_line } , "::end" , NL ,
               "::schema" , NL , { kv_line } , "::end" , NL ,
               "::index"  , NL , { "#" , TEXT_LINE , NL | index_row } , "::end" , NL ;
kv_line      = id , ":" , WS , TEXT_LINE , NL ;
index_row    = id , WS , "|" , WS , node_type , WS , "|" , WS , scope , WS , "|" , WS , TEXT_LINE , NL ;
node_block   = "::node" , WS , "[" , attr , { WS , attr } , "]" , NL ,
               "::dense" , ( WS , TEXT_LINE , NL | NL , TEXT_LINE , NL , "::end" , NL ) ,
               { edge_decl | "::edges" , NL , { edge_decl } , "::end" , NL | fenced_code | colon_block } ,
               "::end" , NL ;
attr         = id , "=" , ( STRING | [ "-" ] , NUMBER | "true" | "false" | id_list ) ;
id_list      = id , { "," , id } ;
node_type    = "figure" | "lattice" | "alphabet" | "frame" | "field" | "embedding" | "invariant" | "meta" ;
scope        = "all" | "clinician" | "orchestrator" | "worker" | id ;
fenced_code  = "```" , ( "hex" | "tla" | "lean" | "mermaid" | id ) , NL , VERBATIM , "```" , NL ;
colon_block  = "::" , id , [ WS , "[" , attr , { WS , attr } , "]" ] , NL , VERBATIM , "::end" , NL ;
edge_decl    = ( "->" | "<-" | "<>" ) , "[:" , edge_type , { WS , attr } , "]" ,
               WS , target , { "," , [ WS ] , target } , NL ;
edge_type    = "contains" | "requires" | "precedes" | "links" | "coarsens" | "supersedes" ;
target       = [ URI_PATH , "#" ] , id ;
id           = ( LETTER | "_" ) , { LETTER | DIGIT | "_" | "-" } ;
field_stream = { record } ;
record       = code_le , word_raw ;
```

## 4. Temporal & Transport Logic (`tla`)

```tla
---- MODULE FigKernel ----
EXTENDS Naturals, Integers, Sequences, FiniteSets

Pow2(n) ==
  LET F[i \in 0..n] == IF i = 0 THEN 1 ELSE 2 * F[i - 1]
  IN F[n]

CodeBytes(d, b) == (d * b + 7) \div 8

Coarsen(z, d, k) == z \div Pow2(k * d)

StrictlyInc(s) ==
  \A i \in 1..(Len(s) - 1) : s[i] < s[i + 1]

CanonicalCodes(s, d, b) ==
  /\ StrictlyInc(s)
  /\ \A i \in 1..Len(s) : 0 <= s[i] /\ s[i] < Pow2(d * b)

CidMatch(fetched, declared) == fetched = declared

NormDisposal(d) == IF d = 0 THEN 1 ELSE d

DrawCanvas(canvas, codes, words, transp) ==
  [z \in DOMAIN canvas |->
    IF z \in codes /\ words[z] /= transp THEN words[z] ELSE canvas[z]]

DisposeCanvas(before, visible, domain, drawn, disp, bg) ==
  CASE NormDisposal(disp) = 1 -> visible
    [] NormDisposal(disp) = 2 ->
         [z \in DOMAIN visible |-> IF z \in domain THEN bg ELSE visible[z]]
    [] NormDisposal(disp) = 3 ->
         [z \in DOMAIN visible |-> IF z \in drawn THEN before[z] ELSE visible[z]]

AutonomousGate(userInput) == userInput = 0
====
```

## 5. Canonical Instantiation (`SYNTH-1937-FIG`)

`SYNTH-1937-FIG (1937-11-14): 2D 4×4 lattice at q = 1 Å (10^-10 m) over GF(2^8) (poly 0x11b), two GIF frames (10 cs keep, 5 cs transparent-hole + restore-background), headerless raw fields pinned by CIDv1 SHA-256.`

````markdown
::meta
title: SYNTH-1937-FIG Multi-Scale Coordinate Telemetry
version: 1.0.0
fig-id: SYNTH-1937-FIG
date: 1937-11-14
scope: all
::end

::schema
node-types: [figure, lattice, alphabet, frame, field, embedding, invariant, meta]
edge-types: [contains, requires, precedes, links, coarsens, supersedes]
scope-levels: [all, clinician, orchestrator, worker]
field-modulus: 0x11b
cid-codec: raw-sha256-cidv1
::end

::index
# id | type | scope | keywords
fig-root | figure | all | loop=1, background=00, centisecond, cidv1
lat-planar | lattice | all | dims=2, bits=2, q=1e-10m, z-order
alph-gf256 | alphabet | all | gf256, poly=0x11b, value-bytes=1
frm-0 | frame | all | seq=0, delay=10, disposal=1, keep
fld-0 | field | all | headerless, 4-records, gf-sum, gf-mul
frm-1 | frame | all | seq=1, delay=5, disposal=2, transparent=00
fld-1 | field | all | headerless, 2-records, hole, sparse
emb-dicom | embedding | all | ipp, iop, spacing, rescale
emb-cif | embedding | all | fract-mod1, cell-metric, angstrom
inv-cid-gate | invariant | all | canonical-codes, cid-match, no-user-input
__view | meta | all | mermaid, stateDiagram-v2
::end

::node [id=fig-root type=figure scope=all loop=1 background="00" time-num=1 time-den=100 var=phase init=init]
::dense SYNTH-1937-FIG|4x4|q=1A|gf256|frames=2|loop=1
->[:contains] lat-planar, alph-gf256, frm-0, frm-1, emb-dicom, emb-cif, inv-cid-gate
->[:supersedes] dicom-pixel-coordinate, cif-atom-site-coordinate
```tla
PhaseDom == {"init", "fr0_held", "fr1_visible", "rested"}
Done == phase = "rested" /\ UNCHANGED vars
```
::end

::node [id=lat-planar type=lattice scope=all dims=2 bits=2 axes=x,y quantum-num=1 quantum-den=10000000000]
::dense dims=2|bits=2|code_bytes=1|q=1e-10m|origin=0,0
```tla
Cells == 0..15
```
::end

::node [id=alph-gf256 type=alphabet scope=all poly=283 value-bytes=1]
::dense gf(2^8)|poly=0x11b|0x57+0x83=0xd4|0x57*0x83=0xc1
```tla
Bg == 0
```
::end

::node [id=frm-0 type=frame scope=all seq=0 delay=10 disposal=1 user-input=0]
::dense seq=0|delay=10cs|disposal=keep|user_input=0
->[:requires] lat-planar, alph-gf256
->[:contains] fld-0
->[:precedes] frm-1
```tla
StepFrame0 ==
  /\ phase = "init" /\ CidMatch(cidOk, TRUE) /\ AutonomousGate(userInput)
  /\ canvas' = DrawCanvas(canvas, F0Codes, F0Words, -1)
  /\ phase' = "fr0_held" /\ UNCHANGED <<cidOk, userInput>>
```
::end

::node [id=fld-0 type=field scope=all records=4 cid="bafkreicefv34tw7gjhd5hikg2tuvg62l3kuzjzpthraj2fnizxgu4gtdcq"]
::dense (0,0)->57|(1,0)->83|(0,1)->d4|(1,1)->c1|token_id=0x442d77c9dbe649c7d3a146d4e9537b4bdaa994e5f33c409d15a8cdcd4e1a6314
```hex
0057018302d403c1
```
::end

::node [id=frm-1 type=frame scope=all seq=1 delay=5 disposal=2 user-input=0 transparent="00"]
::dense seq=1|delay=5cs|disposal=background|transparent=00
->[:requires] lat-planar, alph-gf256
->[:contains] fld-1
::warning
Code 3 stores transparent word 00: cell (1,1) keeps c1 while frm-1 is visible, then disposal=2 restores both stored codes {3,12} to background 00 while unstored code 2 keeps d4.
::end
```tla
StepFrame1 ==
  /\ phase = "fr0_held" /\ CidMatch(cidOk, TRUE) /\ AutonomousGate(userInput)
  /\ canvas' = DrawCanvas(canvas, F1Codes, F1Words, 0)
  /\ phase' = "fr1_visible" /\ UNCHANGED <<cidOk, userInput>>
DisposeFrame1 ==
  /\ phase = "fr1_visible"
  /\ canvas' = DisposeCanvas(canvas, canvas, F1Codes, {12}, 2, Bg)
  /\ phase' = "rested" /\ UNCHANGED <<cidOk, userInput>>
```
::end

::node [id=fld-1 type=field scope=all records=2 cid="bafkreicyel4mobw7ofnc7w4awt74xldi2md5cac6l64heqq3u73s2rad4y"]
::dense (1,1)->00|(2,2)->11|token_id=0x5822f8c706df715a2fdb80b4ffcbac68d307d1005e5fb872421ba7f72d4403e6
```hex
03000c11
```
::end

::node [id=emb-dicom type=embedding scope=all kind=dicom]
::dense IPP=centre(0,0,0)|IOP=(col_cos,row_cos)|cell=roundHalfAway((P_m-origin)/q)|rescale=slope*leInt(w)+intercept
->[:requires] lat-planar, alph-gf256
->[:supersedes] dicom-pixel-coordinate
::end

::node [id=emb-cif type=embedding scope=all kind=cif]
::dense fract_mod_1->metric_A*1e-10->roundHalfAway((X_m-origin)/q)|no-symmetry-orbit
->[:requires] lat-planar, alph-gf256
->[:supersedes] cif-atom-site-coordinate
::end

::node [id=inv-cid-gate type=invariant scope=all kind=safety]
::dense canonical_codes/\cid_match/\user_input=0
->[:requires] fld-0, fld-1, frm-0, frm-1
```tla
InvCidAndGate ==
  /\ CanonicalCodes(<<0, 1, 2, 3>>, 2, 2)
  /\ CanonicalCodes(<<3, 12>>, 2, 2)
  /\ (phase = "fr1_visible" => canvas[3] = 193 /\ canvas[12] = 17)
  /\ (phase = "rested" => canvas[3] = 0 /\ canvas[12] = 0 /\ canvas[2] = 212)
```
::end

::node [id=__view type=meta scope=all]
::dense mermaid|stateDiagram-v2
```mermaid
stateDiagram-v2
  direction LR
  [*] --> Init: bg=00, loop=1
  Init --> Fr0Held: frm-0 [cidOk, delay=10cs, keep {0:57, 1:83, 2:d4, 3:c1}]
  Fr0Held --> Fr1Visible: frm-1 [transp=00 => 3:c1, 12:11, delay=5cs]
  Fr1Visible --> Rested: disposal=2 [wipe {3,12} -> 00; keep {0:57, 1:83, 2:d4}]
  Init --> Halted: cid mismatch or user-input=1 (breaks InvCidAndGate)
  classDef breaks fill:#ffecee,stroke:#e72e3b,stroke-width:2.5px,color:#cd2834
  class Halted breaks
```
::end
````

```mermaid
stateDiagram-v2
  direction LR
  [*] --> Init: bg=00, loop=1
  Init --> Fr0Held: frm-0 [cidOk, delay=10cs, keep {0:57, 1:83, 2:d4, 3:c1}]
  Fr0Held --> Fr1Visible: frm-1 [transp=00 => 3:c1, 12:11, delay=5cs]
  Fr1Visible --> Rested: disposal=2 [wipe {3,12} -> 00; keep {0:57, 1:83, 2:d4}]
  Init --> Halted: cid mismatch or user-input=1 (breaks InvCidAndGate)
  classDef breaks fill:#ffecee,stroke:#e72e3b,stroke-width:2.5px,color:#cd2834
  class Halted breaks
```

## 6. Falsification

Each row is a challenge. The failure signal is what a conforming read returns.

| Challenge | Failure signal |
|---|---|
| Swap axis bit-order on `(3, 2)` at `b = 2` | `mortonEncode 2 [3, 2]` is `13`, not `11`; `(2..3, 2..3)` is the contiguous range `12..15`. |
| Coarsen `(3, 2)` at `b = 4` without shifting `k * d` bits | `coarsen z 2 1` equals `mortonEncode 3 [1, 1]`; fine cell `[x*q, (x+1)*q)` stays inside coarse cell. |
| Separate molecular and anatomical coordinate spaces | At `q = 1 pm`, `b = 48` (`codeBytes 3 48 = 18`), `154` and `1800000000000` round-trip in one lattice and coarsen at `k = 30` to `[1676, 0, 0]`. |
| Sort `code_bytes = 2` records `(1, aa)` and `(256, bb)` lexicographically by raw bytes | Byte-sorted `0001bb0100aa` fails `strictlyIncreasing`; canonical stream is `0100aa0001bb`. |
| Prepend `GIF89a` or `DICM`, duplicate a Morton code, or append a trailing byte | `serializeField` returns `none`; a sparse code is searched by integer key, never multiplied as a byte offset. |
| Compute CIDv1 raw SHA-256 of `""`, `"hello"`, `fld-0`, `fld-1`, or the `SYNTH-1937-FIG` document | Empty is `bafkreihdwdcefgh4dqkjv67uzcmw7ojee6xedzdetojuzjevtenxquvyku`; `"hello"` is `bafkreibm6jg3ux5qumhcn2b3flc3tyu6dmlb4xa7u5bf44yegnrjhc4yeq`; `fld-0` (`0057018302d403c1`) is `bafkreicefv34tw7gjhd5hikg2tuvg62l3kuzjzpthraj2fnizxgu4gtdcq`; `fld-1` (`03000c11`) is `bafkreicyel4mobw7ofnc7w4awt74xldi2md5cac6l64heqq3u73s2rad4y`; `SYNTH-1937-FIG` is `bafkreie2a6pryki5em24exzahodr2l6fgxt5fvo6lyfqaibuvnf27tvlay`. |
| Equate an RFC 6962 chunk Merkle root with the field CID | `fld-0` Merkle root is `c7cd1b18f6821bacf356bedc6d74e3390854da7a2a8fbec1f155d80140f0056f` at 1-record chunks and `e0c8ed099176aed1b9bc45763e82eef4605d212341b5e7020851e324c4c20d56` at 2-record chunks; `cid` stays invariant at `sha256` `442d77c9dbe649c7d3a146d4e9537b4bdaa994e5f33c409d15a8cdcd4e1a6314`. |
| Replace `GF(2^8)` addition or multiplication with mod-256 integer arithmetic | `gfAdd 0x57 0x83 = 0xd4` (`/= 0xda`), `gfMul 0x57 0x83 = 0xc1`, `gfMul 0x57 0xbf = 0x01`; DICOM rescale `rescale [0x90, 0x01] 1 (-1024) = -624` is a separate affine reading. |
| Confuse GIF disposal `2` with disposal `3` on transparent hole `(1, 1) -> 00` | Disposal `2` wipes `{3, 12}` to `00` and preserves code `2 -> d4`; disposal `3` restores only overwritten `{12}` and preserves `3 -> c1`; `bounds = [0, 1]` wipes unstored code `1`. |
| Set `user-input = 1`, omit a play bound on `loop = 0`, or skip `::warning` | `stepFrame` returns `none` on `userInput = true`; unbounded loops require an explicit caller cap; skipping `::warning` invalidates the read. |
| Embed `.dicom` negative `IPP = -10 mm` at `origin = 0`, or quantize `1 mm` spacing at `q = 10 mm` | `dicomVoxel` returns `none` until `origin = -10 mm` (`-> [0, 0, 0]`); colliding cells fail `injectiveCoords`. |
| Embed `.cif` fractional `-1/5` or `1`, or triclinic `CIF-G60` `(a=b=c=2 Å, angles 90,90,60, site (0, 1/2, 0), q=1 pm)` | `-1/5` wraps to `4/5` (`-> 8` on `10 Å`), `1` wraps to `0`, and `CIF-G60` lands at `(50, 87, 0)` since `50 * sqrt(3) pm = sqrt(7500)` is closer to `87` than `86`. |
| Include patient identifiers in `.fig` nodes or treat a passing check as clinical clearance | Any PII attribute is non-conformant; `DISCLAIMER` governs every read. |

`zorder_laws`, `gf256_laws`, `stream_and_gif_laws`, and `embed_laws` are the discharged checks.

## References

1. G. M. Morton, *A Computer Oriented Geodetic Data Base and a New Technique in File Sequencing*, IBM Ltd., 1966.
2. CompuServe Incorporated, [*Graphics Interchange Format Version 89a*](https://www.w3.org/Graphics/GIF/spec-gif89a.txt), §23 Graphic Control Extension, 1990.
3. Netscape Communications, *GIF89a `NETSCAPE2.0` Application Extension* (loop iteration block), 1995.
4. National Institute of Standards and Technology, [*FIPS PUB 197: Advanced Encryption Standard (AES)*](https://doi.org/10.6028/NIST.FIPS.197-upd1), §4.2 `GF(2^8)` arithmetic (`0x11b`), 2001.
5. NEMA, [*DICOM PS3.3: Information Object Definitions*](https://dicom.nema.org/medical/dicom/current/output/html/part03.html), C.7.6.2.1.1 Image Position and Image Orientation (Patient).
6. S. R. Hall, F. H. Allen, and I. D. Brown, ["The Crystallographic Information File (CIF): a new standard archive file for crystallography"](https://doi.org/10.1107/S010876739101067X), *Acta Crystallographica Section A* 47 (1991) 655–685.
7. Protocol Labs, [*CID (Content Identifier) Specification*](https://specs.ipfs.tech/cid/), CIDv1 multicodec `raw` (`0x55`) and multihash `sha2-256` (`0x12`).
8. B. Laurie, A. Langley, and E. Kasper, [*RFC 6962: Certificate Transparency*](https://www.rfc-editor.org/rfc/rfc6962), §2.1 Merkle Hash Trees, 2013.
9. W. Entriken, D. Shirley, J. Evans, and N. Sachs, [*EIP-721: Non-Fungible Token Standard*](https://eips.ethereum.org/EIPS/eip-721), `tokenId` and `tokenURI`, 2018.
10. M. Dubey, [*ObjectGraph* (`arXiv:2604.27820v1`)](https://arxiv.org/html/2604.27820v1), 2026.

## License

MIT — see [LICENSE](LICENSE).
