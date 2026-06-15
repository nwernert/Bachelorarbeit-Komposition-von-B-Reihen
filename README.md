# Bachelorarbeit-Komposition-von-B-Reihen
Dieses Repository enthält code, den ich für meine Bachelorarbeit Komposition von B-Reihen im Sommersemester 2026 an der Johannes Gutenberg Universität Mainz geschrieben habe. Die verwendete Julia Version ist 1.10.5.
Der code soll die Projekte RootedTrees.jl und BSeries.jl erweitern und enthält

-ColoredSplittingIterator, ein Analogon zu SplittingIterator aus RootedTrees.jl für farbige Wurzelbäume.

-compose(b, a, t::BicoloredRootedTree), was es ermöglicht B-Reihen für partitionierte ODEs der Form u'(t) = f(u(t))+g(u(t)) zu komponieren.

-invert(a :: TruncatedBSeries), um B-Reihen basierend auf zwei- und einfarbigen Wurzelbäumen zu invertieren.

-compose(b, a, factor_b, factor_a), eine Verallgemeinerung von compose(b, a; normalize_stepsize = false) für möglicherweise unterschiedliche Schrittweiten.

Um die Tests auszuführen,
