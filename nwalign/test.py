import nwalign as nw
print nw

def test_nw():

    r = nw.global_align("CEOLECANTH", "PELICAN")
    assert r ==  ('CEOLECANTH', 'PE-LICAN--'), r

def test_pam():
    r = nw.global_align("CEELECANTH", "PELICAN", matrix='PAM250')
    assert r ==  ('CEELECANTH', '-PELICA--N')


def test_gap_init():
    s0 = "MIPMDNYCVPSTSTTGLVFSATSSMNASSGFHLTVNSPTSVTGLKHEASLAVDWSVEEQYILEKGLSKFKDEPQVTKYVKIAATLPDKSVRDVAMRCKWMTQKRRKGEEHSTGTKVSYRKVVDLPPKLNMFSTEPQQNATYAMNHMCQSARMPFEGLSDAVMERLRQNAQAFSQISSNLSVCKPQDNVSLFYMARNNISAILNDMKEMPGIISRMPPLPVSINNDLASSLVTSATQPRSYTIPSSIYLKQEPRN"
    s1 = "MIPMDNYCVPSTSTTGLVFSATSSMNASSGFHLTVNSPTSVTGLKHEASLAVDWSVEEQYILEKGLSKFKDEPQVTKYVKIAATLPDKSVRDVAMRCKWMTQKRRKGEEHSTGTKVSYRKVVDLPPKLNMFSTEPQQNATYAMNHMCQSARMPFEGLSDAVMERLRQNAQAFSQISSNLSVCKHEGDAWDHQPDAASACLN"

    a0 = "MIPMDNYCVPSTSTTGLVFSATSSMNASSGFHLTVNSPTSVTGLKHEASLAVDWSVEEQYILEKGLSKFKDEPQVTKYVKIAATLPDKSVRDVAMRCKWMTQKRRKGEEHSTGTKVSYRKVVDLPPKLNMFSTEPQQNATYAMNHMCQSARMPFEGLSDAVMERLRQNAQAFSQISSNLSVCKPQDNVSLFYMARNNISAILNDMKEMPGIISRMPPLPVSINNDLASSLVTSATQPRSYTIPSSIYLKQEPRN"
    a1 = "MIPMDNYCVPSTSTTGLVFSATSSMNASSGFHLTVNSPTSVTGLKHEASLAVDWSVEEQYILEKGLSKFKDEPQVTKYVKIAATLPDKSVRDVAMRCKWMTQKRRKGEEHSTGTKVSYRKVVDLPPKLNMFSTEPQQNATYAMNHMCQSARMPFEGLSDAVMERLRQNAQAFSQISSNLSVC------------------------KHEG--DAWDHQP-----DAASACL----------------------N"

    r = nw.global_align(s0, s1, gap=-2, gap_init=-10, matrix='BLOSUM62')
    assert r[0] == a0
    assert r[1] == a1


if __name__ == "__main__":
    test_gap_init()
