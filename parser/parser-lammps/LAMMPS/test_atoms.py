# Copyright 2015-2018 Massimo Riello, Adam Fekete, Fawzi Mohamed, Ankit Kariryaa
# 
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

import unittest

class TestAtomStyles(unittest.TestCase):

    def test_style_atomic(self):
        style_atom = atom_style('atomic')
        self.assertIsInstance(style_atom, atom_style)

    def test_style_bond(self):
        style_bond = atom_style('bond')
        self.assertIsInstance(style_bond, atom_style)

    def test_style_full(self):
        style_full = atom_style('full')
        self.assertIsInstance(style_full, atom_style)

    def test_style_body(self):
        style_body = atom_style('body', 'nparticle', 2, 10)
        self.assertIsInstance(style_body, atom_style)

    def test_style_hybrid(self):
        style_hybrid = atom_style('hybrid', 'charge', 'bond')
        self.assertIsInstance(style_hybrid, atom_style)

    def test_style_hybrid(self):
        style_hybrid = atom_style('hybrid', 'charge', 'body', 'nparticle', 2, 5)
        self.assertIsInstance(style_hybrid, atom_style)

    def test_style_template(self):
        style_template = atom_style('template', 'myMols')
        self.assertIsInstance(style_template, atom_style)

