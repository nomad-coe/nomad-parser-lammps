package eu.nomad_lab.parsers

import org.specs2.mutable.Specification

object LammpsParserSpec extends Specification {
  "LammpsParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(LammpsParser, "parsers/lammps/test/examples/methane/input.methane_nvt_thermo_style_custom", "json-events") must_== ParseResult.ParseSuccess
    }
  }
}