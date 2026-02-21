#!/usr/bin/env python
__copyright__ = (
    "Copyright (©) 2020-2025 EPFL (Ecole Polytechnique Fédérale de Lausanne) "
    "Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)"
)
__license__ = "LGPLv3"

import math

import akantu as aka
import numpy as np
import pytest


class TestParameterRegistry:
    """Tests for the ParameterRegistry interface"""

    def test_register_real_without_default(self):
        """Test registering a Real parameter without default value"""
        parsable = aka.Parsable(aka.ParserType.material, "test_material")

        # Register a parameter without default
        parsable.registerParamReal("E", aka._pat_modifiable, "Young's modulus")

        # Set and get the value
        parsable.set("E", 210000.0)
        value = parsable.get("E")
        assert value == pytest.approx(210000.0)

    def test_register_real_with_default(self):
        """Test registering a Real parameter with default value"""
        parsable = aka.Parsable(aka.ParserType.material, "test_material")

        # Register a parameter with default
        parsable.registerParamReal(100.0, "E", aka._pat_modifiable, "Density")

        # Get default value
        value = parsable.get("E")
        assert value == pytest.approx(100.0)


class TestParser:
    """Tests for the Parser class"""

    def test_parse_input(self):
        """Test parsing an existing input file"""
        aka.parseInput("input_file.dat")

        parser = aka.getStaticParser()

        param = parser.getParameter("toto")
        assert param.getValue() == "2*pi + max(2, general)"

        v = param.asReal()
        assert v == pytest.approx(2 * math.pi + max(2, 50))

        mat = parser.getParameter("mat").asMatrix()
        ref = np.array([[1, 25, 5, v], [0, 10, 50, 13]])
        assert np.linalg.norm(mat - ref, np.inf) == pytest.approx(0)


class TestParsable:
    """Tests for the Parser class"""

    def test_parsable(self):
        """Test registering a Real parameter without default value"""
        parsable = aka.Parsable(aka.ParserType.material, "test_material")

        # Register a parameter without default
        parsable.registerParamReal("E", aka._pat_parsmod, "Young's modulus")
        parsable.registerParamReal(
            "X135", aka._pat_parsable, "Random parameter"
        )

        parser = aka.getStaticParser()

        sub_section = parser.getSubSections(aka.ParserType.material)

        E = sub_section.getParameter("E").asReal()
        assert E == pytest.approx(1)

        parsable.parseSection(sub_section)

        E = parsable.get("E")
        assert E == pytest.approx(1)

        with pytest.raises(aka.Exception):
            parsable.get("X135")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
