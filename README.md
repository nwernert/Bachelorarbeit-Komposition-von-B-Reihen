# Bachelorarbeit-Komposition-von-B-Reihen
Dieses Repository enthält code, den ich für meine Bachelorarbeit Komposition von B-Reihen im Sommersemester 2026 an der Johannes Gutenberg Universität Mainz geschrieben habe. Die verwendete Julia Version ist 1.10.5.
Der code soll die Projekte RootedTrees.jl und BSeries.jl erweitern und enthält

-ColoredSplittingIterator, ein Analogon zu SplittingIterator aus RootedTrees.jl für farbige Wurzelbäume,

-compose(b, a, t::BicoloredRootedTree), was es ermöglicht B-Reihen für partitionierte ODEs der Form u'(t) = f(u(t))+g(u(t)) zu komponieren,

-invert(a :: TruncatedBSeries), um B-Reihen basierend auf zwei- und einfarbigen Wurzelbäumen zu invertieren,

-compose(b, a, factor_b, factor_a), eine Verallgemeinerung von compose(b, a; normalize_stepsize = false) für möglicherweise unterschiedliche Schrittweiten.

Die obigen Funktionen und Tests für sie sind in der Datei invert_and_compose.jl enthalten. Um die Tests in dem Julia REPL auszuführen, führt man folgende Schritte aus:
1. Mit cd() navigiert man sich in den Ordner, wo Manifest.toml, Project.toml und invert_and_compose.jl enthalten sind.
2. Mit Pkg.activate(".") aktiviert man das Projekt mit der richtigen Umgebung.
3. Mit include("invert_and_compose.jl") führt man die Tests aus.
