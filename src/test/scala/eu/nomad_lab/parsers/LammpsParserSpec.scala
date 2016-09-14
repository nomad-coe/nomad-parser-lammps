package eu.nomad_lab.parsers

import org.specs2.mutable.Specification

object LammpsParserSpec extends Specification {
  "LammpsParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(LammpsParser, "parsers/lammps/test/examples/methane/log.64xmethane-nvt-thermo_style_custom", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(LammpsParser, "parsers/lammps/test/examples/methane/log.64xmethane-nvt-thermo_style_custom", "json") must_== ParseResult.ParseSuccess
    }
  }
}

object LammpsParserSpec1 extends Specification {
  "LammpsParserTest1" >> {
    "test with json-events 1" >> {
      ParserRun.parse(LammpsParser, "parsers/lammps/test/examples/methane/log.64xmethane-nvt-thermo_style_one", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json 1" >> {
      ParserRun.parse(LammpsParser, "parsers/lammps/test/examples/methane/log.64xmethane-nvt-thermo_style_one", "json") must_== ParseResult.ParseSuccess
    }
  }
}

object LammpsParserSpec2 extends Specification {
  "LammpsParserTest2" >> {
    "test with json-events 2" >> {
      ParserRun.parse(LammpsParser, "parsers/lammps/test/examples/methane/log.64xmethane-nvt-thermo_style_multi", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json 2" >> {
      ParserRun.parse(LammpsParser, "parsers/lammps/test/examples/methane/log.64xmethane-nvt-thermo_style_multi", "json") must_== ParseResult.ParseSuccess
    }
  }
}

// object LammpsParserSpec3 extends Specification {
//   "LammpsParserTest3" >> {
//     "test with json-events 3" >> {
//       ParserRun.parse(LammpsParser, "parsers/lammps/test/examples/1_methyl_naphthalene/log.1_methyl_naphthalene", "json-events") must_== ParseResult.ParseSuccess
//     }
//     "test with json 3" >> {
//       ParserRun.parse(LammpsParser, "parsers/lammps/test/examples/1_methyl_naphthalene/log.1_methyl_naphthalene", "json") must_== ParseResult.ParseSuccess
//     }
//   }
// }
//
object LammpsParserSpec4 extends Specification {
  "LammpsParserTest4" >> {
    "test with json-events 4" >> {
      ParserRun.parse(LammpsParser, "parsers/lammps/test/examples/hexane_cyclohexane/log.hexane_cyclohexane_nvt", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json 4" >> {
      ParserRun.parse(LammpsParser, "parsers/lammps/test/examples/hexane_cyclohexane/log.hexane_cyclohexane_nvt", "json") must_== ParseResult.ParseSuccess
    }
  }
}
//
//object LammpsParserSpec5 extends Specification {
//  "LammpsParserTest5" >> {
//    "test with json-events 5" >> {
//      ParserRun.parse(LammpsParser, "parsers/lammps/test/examples/methane_MD_traj/log.methane_nvt_traj_atom_thermo_style_custom", "json-events") must_== ParseResult.ParseSuccess
//    }
//    "test with json 5" >> {
//      ParserRun.parse(LammpsParser, "parsers/lammps/test/examples/methane_MD_traj/log.methane_nvt_traj_atom_thermo_style_custom", "json") must_== ParseResult.ParseSuccess
//    }
//  }
//}
//
//object LammpsParserSpec6 extends Specification {
//  "LammpsParserTest6" >> {
//    "test with json-events 6" >> {
//      ParserRun.parse(LammpsParser, "parsers/lammps/test/examples/methane_MD_traj/log.methane_nvt_traj_atom_unscaled_image_flags_thermo_style_custom", "json-events") must_== ParseResult.ParseSuccess
//    }
//    "test with json 6" >> {
//      ParserRun.parse(LammpsParser, "parsers/lammps/test/examples/methane_MD_traj/log.methane_nvt_traj_atom_unscaled_image_flags_thermo_style_custom", "json") must_== ParseResult.ParseSuccess
//    }
//  }
//}
//
//object LammpsParserSpec7 extends Specification {
//  "LammpsParserTest7" >> {
//    "test with json-events 7" >> {
//      ParserRun.parse(LammpsParser, "parsers/lammps/test/examples/methane_MD_traj/log.methane_nvt_traj_lammpstrj_thermo_style_custom", "json-events") must_== ParseResult.ParseSuccess
//    }
//    "test with json 7" >> {
//      ParserRun.parse(LammpsParser, "parsers/lammps/test/examples/methane_MD_traj/log.methane_nvt_traj_lammpstrj_thermo_style_custom", "json") must_== ParseResult.ParseSuccess
//    }
//  }
//}
//
//object LammpsParserSpec8 extends Specification {
//  "LammpsParserTest8" >> {
//    "test with json-events 8" >> {
//      ParserRun.parse(LammpsParser, "parsers/lammps/test/examples/methane_MD_traj/log.methane_nvt_traj_xyz_thermo_style_custom", "json-events") must_== ParseResult.ParseSuccess
//    }
//    "test with json 8" >> {
//      ParserRun.parse(LammpsParser, "parsers/lammps/test/examples/methane_MD_traj/log.methane_nvt_traj_xyz_thermo_style_custom", "json") must_== ParseResult.ParseSuccess
//    }
//  }
//}