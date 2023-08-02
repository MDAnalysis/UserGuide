
from doc.source.scripts.gen_unit_tables import write_unit_table
import tempfile

def test_write_unit_table(snapshot):
    with tempfile.NamedTemporaryFile() as f:
        tables = write_unit_table(f.name)
        assert tables == snapshot
