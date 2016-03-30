package eu.nomad_lab.parsers

import eu.{ nomad_lab => lab }
import eu.nomad_lab.DefaultPythonInterpreter
import org.{ json4s => jn }
import scala.collection.breakOut

object LammpsParser extends SimpleExternalParserGenerator(
  name = "LammpsParser",
  parserInfo = jn.JObject(
    ("name" -> jn.JString("LammpsParser")) ::
      ("parserId" -> jn.JString("LammpsParser" + lab.LammpsVersionInfo.version)) ::
      ("versionInfo" -> jn.JObject(
        ("nomadCoreVersion" -> jn.JString(lab.NomadCoreVersionInfo.version)) ::
          (lab.LammpsVersionInfo.toMap.map {
            case (key, value) =>
              (key -> jn.JString(value.toString))
          }(breakOut): List[(String, jn.JString)])
      )) :: Nil
  ),
  mainFileTypes = Seq("text/.*"),
  mainFileRe = """# input script for topotools""".r,
  cmd = Seq(DefaultPythonInterpreter.python2Exe(), "${envDir}/parsers/lammps/parser/parser-lammps/LAMMPSParserControl.py",
    "${mainFilePath}"),
  resList = Seq(
    "parser-lammps/LAMMPSParserControl.py",
    "parser-lammps/LAMMPSParserData.py",
    "parser-lammps/LAMMPSParserInput.py",
    "parser-lammps/LAMMPSParserLog.py",
    "parser-lammps/setup_paths.py",
    "nomad_meta_info/common.nomadmetainfo.json",
    "nomad_meta_info/meta_types.nomadmetainfo.json",
    "nomad_meta_info/lammps.nomadmetainfo.json"
  ) ++ DefaultPythonInterpreter.commonFiles(),
  dirMap = Map(
    "parser-lammps" -> "parsers/lammps/parser/parser-lammps",
    "nomad_meta_info" -> "nomad-meta-info/meta_info/nomad_meta_info"
  ) ++ DefaultPythonInterpreter.commonDirMapping()
)
