import pytest

from rxnchirality import ReactionParseError, Status, analyze_reaction

from helpers import load_gold, oracle_cip

GOLD = load_gold()


@pytest.mark.parametrize("name", sorted(GOLD))
def test_gold(name):
    case = GOLD[name]
    res = analyze_reaction(case["rxn"])
    got = {
        c.map_num: (c.status.value, c.reactant_cip, c.product_cip, c.product_rcip)
        for c in res.centers
    }
    assert got == case["centers"]


@pytest.mark.parametrize("name", sorted(GOLD))
def test_gold_labels_match_rdkit(name):
    reactants, _, products = GOLD[name]["rxn"].split(">")
    for m, (_, r_cip, p_cip, _) in GOLD[name]["centers"].items():
        assert r_cip == oracle_cip(reactants, m)
        assert p_cip == oracle_cip(products, m)


@pytest.mark.parametrize("name", sorted(GOLD))
def test_rcip_flips_exactly_on_inversion(name):
    for c in analyze_reaction(GOLD[name]["rxn"]).centers:
        if c.status in (Status.RETAINED, Status.INVERTED) and c.product_rcip:
            assert (c.product_rcip != c.reactant_cip) == (c.status is Status.INVERTED)


def test_counts():
    res = analyze_reaction(GOLD["two_centers"]["rxn"])
    assert (res.n_inverted, res.n_retained, res.n_created, res.n_ambiguous) == (1, 1, 0, 0)
    assert res.summary()["inverted"] == 1


def test_slot_inheritance_ligand_sets():
    c = analyze_reaction(GOLD["phosphine_oxidation_retention"]["rxn"]).center(2)
    assert c.leaving == ("LP",) and c.entering == (15,)
    assert c.conserved == (1, 3, 5)


def test_guard_off_accepts_mapped_pseudo_center():
    rxn = GOLD["pseudo_center_guarded"]["rxn"]
    guarded = analyze_reaction(rxn)
    assert guarded.centers == []
    assert {s.reason for s in guarded.skipped} == {"tagged_but_not_stereogenic"}
    c = analyze_reaction(rxn, require_stereogenic=False).center(2)
    assert c.status is Status.INVERTED
    assert c.reactant_cip is None and c.reason == "cip_tie_rules_1_2"


def test_unmapped_center_is_skipped():
    res = analyze_reaction("C[C@H](O)CC>>C[C@H](OC(C)=O)CC")
    assert res.centers == []
    assert [s.reason for s in res.skipped] == ["center_atom_unmapped"] * 2


def test_agents_are_ignored():
    rxn = "[CH3:1][C@H:2]([OH:3])[CH2:4][CH3:5]>[Na+].[OH-]>[CH3:1][C@@H:2]([OH:3])[CH2:4][CH3:5]"
    assert analyze_reaction(rxn).center(2).status is Status.INVERTED


@pytest.mark.parametrize(
    "rxn",
    [
        "CCO>CC",
        "[CH3:1][CH3:1]>>[CH3:1]",
        "[CH3:1][OH:2]>>[CH3:1][NH2:2]",
        "C1CC>>CC",
    ],
)
def test_parse_errors(rxn):
    with pytest.raises(ReactionParseError):
        analyze_reaction(rxn)
