"""End-to-end checks of extract_geometry.py on small .scad models with known answers.

Run: uv run --with pytest pytest tests (needs openscad on PATH)."""
import json
import pathlib
import shutil
import subprocess

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
TESTS = ROOT / "tests"

pytestmark = pytest.mark.skipif(shutil.which("openscad") is None, reason="openscad not on PATH")


def run_raw(cwd, *args):
    """extract_geometry.py in cwd, where .scad exports land."""
    return subprocess.run(["uv", "run", "-q", str(ROOT / "extract_geometry.py"), *map(str, args)],
                          cwd=cwd, capture_output=True, text=True)


def run(cwd, *args):
    proc = run_raw(cwd, *args)
    assert proc.returncode == 0, proc.stderr
    return json.loads(proc.stdout)


def run_brief(cwd, *args):
    proc = run_raw(cwd, *args, "--brief")
    assert proc.returncode == 0, proc.stderr
    return proc.stdout


def thread(cwd, *args):
    return run(cwd, "thread", TESTS / "thread_bore.scad", "--center", "0,0,0", "--rmax", "8", *args)


def assert_known_thread(t, hand="right"):
    assert t["pitch"] == pytest.approx(3.125, abs=0.01)
    assert t["lead"] == pytest.approx(3.125, abs=0.02)
    assert t["starts"] == 1
    assert t["handedness"] == hand
    assert t["minor_radius"] == pytest.approx(6.75, abs=0.01)
    assert t["major_radius"] == pytest.approx(7.15, abs=0.01)
    assert t["root_flat"] == pytest.approx(0.625, abs=0.03)
    # the groove centre sits half a pitch from a crest centre
    phase = (t["crest_center_at_0deg"] + 3.125 / 2 - 5.44375) % 3.125
    assert min(phase, 3.125 - phase) < 0.05


def test_thread_round_bore(tmp_path):
    assert_known_thread(thread(tmp_path))


def test_thread_left_hand(tmp_path):
    t = thread(tmp_path, "-D", "left=true")
    assert t["handedness"] == "left"
    assert t["pitch"] == pytest.approx(3.125, abs=0.01)


def test_thread_octagonal_bore_with_slot(tmp_path):
    # with many rays, those near the octagon's corners see a constant radius between the
    # root and --rmax; they must not outvote the root (this once gave major 7.18, null pitch)
    assert_known_thread(thread(tmp_path, "-D", "oct=true", "--rmax", "7.3", "--angles", "72"))


def test_thread_plain_bore_is_an_error(tmp_path):
    proc = run_raw(tmp_path, "thread", TESTS / "thread_bore.scad", "--center", "0,0,0", "--rmax", "8",
                   "-D", "threaded=false")
    assert proc.returncode != 0
    assert "no thread" in proc.stderr


def test_ascii_stl(tmp_path):
    stl = tmp_path / "cube.stl"
    (tmp_path / "cube.scad").write_text("cube(10);")
    subprocess.run(["openscad", "--export-format", "asciistl", "-o", stl, tmp_path / "cube.scad"],
                   check=True, capture_output=True)
    assert stl.read_bytes().startswith(b"solid")
    stats = run(tmp_path, stl, "--stats")
    assert stats["volume"] == pytest.approx(1000)
    assert stats["watertight"]


@pytest.mark.parametrize("define, scale", [("turn=true", "-1,-1,1"), ("mirrored=true", "-1,1,1")])
def test_compare_scale_undoes_a_transform(tmp_path, define, scale):
    block = TESTS / "l_block.scad"
    d = run(tmp_path, "compare", block, block, "-Db", define, "--scale", scale, "--slice", "Z=2.5")
    assert d["b_stats"]["volume_sign"] == 1
    assert d["volume_diff"] == pytest.approx(0, abs=1e-3)
    assert all(abs(r["diff"]) < 1e-3 for r in d["area_by_normal"])
    s = d["slices"][0]
    assert not s["unmatched_a"] and not s["unmatched_b"] and not s["loop_diffs"]


def test_compare_without_scale_sees_the_transform(tmp_path):
    block = TESTS / "l_block.scad"
    d = run(tmp_path, "compare", block, block, "-Db", "turn=true", "--slice", "Z=2.5")
    assert d["slices"][0]["unmatched_a"]


def test_brief_faces_slice_and_compare(tmp_path):
    block = TESTS / "l_block.scad"
    text = run_brief(tmp_path, block)
    assert text.startswith("volume 1187.")  # 20x10x5 + 5x10x5 - a r2 hole
    assert "hole centres (12, 5, 0)" in text or "hole centres (12, 5, 5)" in text

    text = run_brief(tmp_path, "slice", block, "--origin", "0,0,2.5", "--normal", "0,0,1")
    circles = [line.split() for line in text.splitlines() if line.strip().startswith("C ")]
    assert len(circles) == 1 and float(circles[0][-1]) == pytest.approx(2, abs=0.01)
    assert any(line.startswith("hole ") for line in text.splitlines())

    text = run_brief(tmp_path, "compare", block, block, "-Db", "turn=true", "--slice", "Z=2.5")
    assert "volume_diff 0" in text
    assert "only_a outer" in text
