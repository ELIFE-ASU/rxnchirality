import json

from rxnchirality.cli import main

from helpers import load_gold

GOLD = load_gold()


def test_tsv_centers(capsys):
    assert main([GOLD["two_centers"]["rxn"]]) == 0
    lines = capsys.readouterr().out.strip().splitlines()
    assert lines[0].split("\t")[:4] == ["id", "map_num", "element", "status"]
    rows = [line.split("\t") for line in lines[1:]]
    assert [(r[1], r[3]) for r in rows] == [("2", "retained"), ("4", "inverted")]


def test_summary_from_file(tmp_path, capsys):
    f = tmp_path / "rxns.smi"
    f.write_text(f"# comment\n{GOLD['sn2']['rxn']}\tsn2\nnot>a>reaction>x\tbad\n")
    assert main(["-i", str(f), "--summary"]) == 1
    lines = capsys.readouterr().out.strip().splitlines()
    header = lines[0].split("\t")
    sn2 = dict(zip(header, lines[1].split("\t")))
    assert sn2["id"] == "sn2" and sn2["inverted"] == "1" and sn2["retained"] == "0"
    assert lines[2].split("\t")[0] == "bad"


def test_json(capsys):
    assert main([GOLD["mitsunobu"]["rxn"], "--format", "json"]) == 0
    data = json.loads(capsys.readouterr().out)
    assert data[0]["summary"]["inverted"] == 1
    assert data[0]["centers"][0]["product_rcip"] == "R"
