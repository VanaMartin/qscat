"""The selected-property run changes only the predeclared sector counts."""

import pytest

from projects.ukrmol_co.fifty_property_export import NEW_SELECTION, OLD_SELECTION, selected_deck


def test_only_the_known_native_selection_changes():
    text = "&input\n ipol=2,\n" + OLD_SELECTION + " ! retained roots\n zlast=.true.,\n/\n"
    result = selected_deck(text)
    assert result.replace(NEW_SELECTION, OLD_SELECTION) == text
    assert NEW_SELECTION in result


@pytest.mark.parametrize("text", ["ntgtl=7,6,6,8,6,8,3,6,", OLD_SELECTION * 2, "ntgtl=5,"])
def test_unknown_or_duplicate_native_selection_is_rejected(text):
    with pytest.raises(ValueError):
        selected_deck(text)
