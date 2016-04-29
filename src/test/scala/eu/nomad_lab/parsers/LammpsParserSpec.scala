package eu.nomad_lab.parsers

import eu.nomad_lab.{ parsers, DefaultPythonInterpreter }
import org.scalacheck.Properties
import org.specs2.mutable.Specification
import org.{ json4s => jn }

object LammpsParserSpec extends Specification {
  “LammpsParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(VaspParser, "test/examples/methane/input.methane_nvt_thermo_style_custom", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(VaspParser, "test/examples/methane/input.methane_nvt_thermo_style_custom", "json") must_== ParseResult.ParseSuccess
    }
  }
}
