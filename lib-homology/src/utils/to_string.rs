use std::fmt::Display;

pub trait ToSeqString<I: Display> {
    fn label(&self) -> String;
    fn indices(&self) -> Vec<I>;
    fn entry_at(&self, i: &I) -> String;
    fn to_seq_string(&self) -> String {
        use yui_core::util::format::table;
        table(self.label(), [""], self.indices(), |_, i| {
            self.entry_at(i)
        })
    }
}

pub trait ToTableString<I: Display> {
    fn labels(&self) -> (String, String);
    fn indices(&self) -> (Vec<I>, Vec<I>);
    fn entry_at(&self, i: &I, j: &I) -> String;
    fn to_table_string(&self) -> String {
        use yui_core::util::format::table;

        let (label0, label1) = self.labels();
        let (ind0, ind1) = self.indices();
        let head = format!("{}\\{}", label1, label0);

        table(head, ind1.into_iter().rev(), ind0, |j, i| {
            self.entry_at(i, j)
        })
    }
}

