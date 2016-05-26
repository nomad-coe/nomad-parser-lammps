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

object LammpsParserSpec extends Specification {
  "LammpsParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(LammpsParser, "parsers/lammps/test/examples/methane/input.methane_nvt_thermo_style_one", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(LammpsParser, "parsers/lammps/test/examples/methane/input.methane_nvt_thermo_style_one", "json") must_== ParseResult.ParseSuccess
    }
  }
}

object LammpsParserSpec extends Specification {
  "LammpsParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(LammpsParser, "parsers/lammps/test/examples/methane/input.methane_nvt_thermo_style_multi", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(LammpsParser, "parsers/lammps/test/examples/methane/input.methane_nvt_thermo_style_multi", "json") must_== ParseResult.ParseSuccess
    }
  }
}

object LammpsParserSpec extends Specification {
  "LammpsParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(LammpsParser, "parsers/lammps/test/examples/1_methyl_naphthalene/input.1_methyl_naphthalene", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(LammpsParser, "parsers/lammps/test/examples/1_methyl_naphthalene/input.1_methyl_naphthalene", "json") must_== ParseResult.ParseSuccess
    }
  }
}