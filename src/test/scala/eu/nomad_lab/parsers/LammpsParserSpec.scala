package eu.nomad_lab.parsers

import org.specs2.mutable.Specification

object LammpsParserSpec extends Specification {
  "LammpsParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(LammpsParser, "parsers/lammps/test/examples/methane/input.methane_nvt_thermo_style_custom", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(LammpsParser, "parsers/lammps/test/examples/methane/input.methane_nvt_thermo_style_custom", "json") must_== ParseResult.ParseSuccess
    }
  }
}

object LammpsParserSpec1 extends Specification {
  "LammpsParserTest1" >> {
    "test with json-events 1" >> {
      ParserRun.parse(LammpsParser, "parsers/lammps/test/examples/methane/input.methane_nvt_thermo_style_one", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json 1" >> {
      ParserRun.parse(LammpsParser, "parsers/lammps/test/examples/methane/input.methane_nvt_thermo_style_one", "json") must_== ParseResult.ParseSuccess
    }
  }
}

object LammpsParserSpec2 extends Specification {
  "LammpsParserTest2" >> {
    "test with json-events 2" >> {
      ParserRun.parse(LammpsParser, "parsers/lammps/test/examples/methane/input.methane_nvt_thermo_style_multi", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json 2" >> {
      ParserRun.parse(LammpsParser, "parsers/lammps/test/examples/methane/input.methane_nvt_thermo_style_multi", "json") must_== ParseResult.ParseSuccess
    }
  }
}

