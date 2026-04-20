"""Unit tests for :mod:`scripts.common.exceptions`.

Verifies the inheritance contract documented in the module docstring:
:class:`InvalidInput` is a ``ValueError`` so every existing
``pytest.raises(ValueError)`` test still matches after the promotion.
"""

import pytest


def test_biko_error_is_exception():
    from scripts.common.exceptions import BikoError

    assert issubclass(BikoError, Exception)


def test_invalid_input_is_value_error_and_biko_error():
    from scripts.common.exceptions import BikoError, InvalidInput

    assert issubclass(InvalidInput, ValueError)
    assert issubclass(InvalidInput, BikoError)


def test_malformed_vcf_is_invalid_input():
    from scripts.common.exceptions import InvalidInput, MalformedVCF

    assert issubclass(MalformedVCF, InvalidInput)
    assert issubclass(MalformedVCF, ValueError)


def test_database_unavailable_is_biko_error_not_value_error():
    from scripts.common.exceptions import BikoError, DatabaseUnavailable

    assert issubclass(DatabaseUnavailable, BikoError)
    assert not issubclass(DatabaseUnavailable, ValueError)


def test_external_api_error_is_biko_error_not_value_error():
    from scripts.common.exceptions import BikoError, ExternalAPIError

    assert issubclass(ExternalAPIError, BikoError)
    assert not issubclass(ExternalAPIError, ValueError)


def test_legacy_value_error_catch_still_matches_invalid_input():
    """The whole point of the multi-inheritance is that old tests keep working."""
    from scripts.common.exceptions import InvalidInput

    with pytest.raises(ValueError):
        raise InvalidInput("demo")


def test_biko_error_catch_matches_any_domain_subclass():
    from scripts.common.exceptions import (
        BikoError,
        DatabaseUnavailable,
        ExternalAPIError,
        InvalidInput,
        MalformedVCF,
    )

    for cls in (InvalidInput, MalformedVCF, DatabaseUnavailable, ExternalAPIError):
        with pytest.raises(BikoError):
            raise cls("demo")


def test_parse_ped_malformed_raises_invalid_input(tmp_path):
    """Regression: PED malformed-row raise site now returns InvalidInput
    but legacy ``except ValueError`` continues to match."""
    from scripts.common.exceptions import InvalidInput
    from scripts.intake.parse_ped import parse_ped

    ped = tmp_path / "bad.ped"
    ped.write_text("FAM1\tCHILD\tFATHER\tMOTHER\t1\t2\nFAM1\tONLY3COLS\tFATHER\n")

    with pytest.raises(InvalidInput):
        parse_ped(str(ped))

    # And legacy ValueError also catches it:
    with pytest.raises(ValueError):
        parse_ped(str(ped))


def test_parse_vcf_strict_ped_raises_invalid_input(tmp_path):
    """Regression: the strict-mode PED resolution failure in parse_vcf
    now raises InvalidInput (ValueError-compatible)."""
    from scripts.common.exceptions import InvalidInput
    from scripts.intake.parse_vcf import parse_vcf

    ped = tmp_path / "orphan.ped"
    ped.write_text("FAM1\tX\tY\tZ\t1\t2\n")
    vcf = tmp_path / "sample.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tFOO\tBAR\tBAZ\n"
        "chr1\t100\t.\tA\tT\t100\tPASS\t.\tGT\t0/0\t0/0\t0/1\n"
    )

    with pytest.raises(InvalidInput):
        parse_vcf(str(vcf), ped_path=str(ped))
