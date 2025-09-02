import csv

class IEDBEpitopeTableHandler:
    def __init__(self, csv_file):
        self.rows = []
        self.header = []
        with open(csv_file, newline='') as f:
            reader = csv.reader(f)
            self.header = next(reader)
            for row in reader:
                self.rows.append(row)

    def get_by_epitope_id(self, epitope_id):
        """
        Returns the row (as a dict) for the given Epitope_id (first column).
        """
        epitope_id_col = self.header.index('Epitope_id')
        for row in self.rows:
            if row[epitope_id_col] == str(epitope_id):
                return dict(zip(self.header, row))
        return None


