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

