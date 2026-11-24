"""The examples: each fills the form, builds, and shows what its name says."""

from dataclasses import replace

import pytest
from gui_fakes import FakeScheduler, FakeView, FakeWorker

from eulerase.gui import presenter as presenter_module
from eulerase.gui.presenter import EXAMPLES, READY, FormValues, Presenter

READS = ">read1\nACGTACGGA\n>read2\nTTGCCATTA"
OVERLAPPING = ">r1\nACGTGCAATG\n>r2\nCAATGCCGTA\n>r3\nCCGTAACGTG"
WITH_ERROR = ">ok\nACGTGCGTACGT\n>err\nACGTTT"

# The forms the eight examples must put in the controls, typed from the spec's table and
# not read from EXAMPLES: a changed text, mode, k or flag is caught here.
# name: (text, input_mode, k_text, both_strands, clip_tips)
EXPECTED = {
    "Lemon text, k=3 (single path)": ("LEMONLEMONADE", "sequences", "3", False, False),
    "DNA with a repeat, k=4": ("ACGTTGCAACGTTGCA", "sequences", "4", False, False),
    "Tandem repeat (self-loop), k=3": ("AAAAAC", "sequences", "3", False, False),
    "Two unrelated reads (2 trails), k=4": (READS, "sequences", "4", False, False),
    "Both strands, k=5": ("ATGGCGTACGTTAGC", "sequences", "5", True, False),
    "A k-mer list that closes into a circuit": ("ACG CGT GTA TAC", "kmers", "3", False, False),
    "Overlapping reads (assembly), k=5": (OVERLAPPING, "sequences", "5", False, False),
    "Reads with one error, tip clipping on, k=4": (WITH_ERROR, "sequences", "4", False, True),
}
NAMES = list(EXPECTED)

# name: a fact about the built result that is independent of how the example is written.
FACTS = {
    NAMES[0]: lambda a: a.headline.startswith("Eulerian path"),
    NAMES[3]: lambda a: len(a.trails) == 2,
    NAMES[5]: lambda a: a.headline.startswith("Eulerian circuit"),
    NAMES[7]: lambda a: (a.clipped_branches, a.clipped_kmers) == (1, 2),
}


def rig():
    view = FakeView()
    view.k_enabled = True  # so a disabled box can only come from the presenter
    worker = FakeWorker(view.calls)
    return view, worker, Presenter(view, worker, FakeScheduler())


def built(name):
    view, worker, presenter = rig()
    presenter.load_example(name)
    worker.finish(0, presenter)
    return view, presenter


def form_of(name):
    text, mode, k_text, both_strands, clip_tips = EXPECTED[name]
    return FormValues(text, mode, k_text, both_strands, False, "", clip_tips)


def test_the_examples_are_the_eight_of_the_spec_in_order():
    assert [e.name for e in EXAMPLES] == NAMES


@pytest.mark.parametrize("name", NAMES)
def test_an_example_builds_and_shows_what_its_name_says(name):
    view, presenter = built(name)
    assert (presenter.phase, view.status_is_error) == (READY, False)
    assert FACTS.get(name, lambda analysis: True)(presenter.analysis)


@pytest.mark.parametrize("name", NAMES)
def test_loading_an_example_fills_the_form_and_builds(name):
    view, worker, presenter = rig()
    assert view.form != form_of(name)  # the fake starts with another form

    presenter.load_example(name)

    assert view.form == form_of(name)
    assert view.k_enabled is (view.form.input_mode == "sequences")
    assert len(worker.jobs) == 1


def test_loading_an_example_carries_its_unique_flag_and_minimum_into_the_form(monkeypatch):
    # No shipped example sets either, so a row that does is put in for this one case.
    row = replace(EXAMPLES[0], unique=True, min_multiplicity=2)
    monkeypatch.setattr(presenter_module, "EXAMPLES", (row,))
    view, _, presenter = rig()

    presenter.load_example(row.name)

    assert view.form == FormValues("LEMONLEMONADE", "sequences", "3", False, True, "2", False)


def test_open_with_none_loads_the_first_example_and_a_form_is_built_as_given():
    given = FormValues("ACG CGT GTA TAC", "kmers", "7", True, True, "2", False)
    assert given != form_of(NAMES[0])

    view, worker, presenter = rig()
    presenter.open_with(given)
    assert view.form == given
    assert view.k_enabled is False

    presenter.open_with(None)
    assert view.form == form_of(NAMES[0])
    assert view.k_enabled is True
    assert len(worker.jobs) == 2


def test_no_example_names_a_parked_feature():
    assert not [e.name for e in EXAMPLES if "binning" in (e.name + e.text).lower()]
