Here’s a list proteins that will give good hits for the spreadsheets that I’ve sent yesterday. The other thing about these ‘good hits’ proteins is that there are regions that would give very bad score/negative hits so we could immediately tell if there’s something not quite right.

Good hits
> Human preproinsulin
MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN

> Human IAPP
MGILKLQVFLIVLSVALNHLKATPIESHQVEKRKCNTATCATQRLANFLVHSSNNFGAILSSTNVGSNTYGKRNAVEVLKREPLNYLPL

>Human IGRP
ALSHTVCGMDKFSITLHRLTWSFLWSVFWLIQISVCISRVFIATHFPHQVILGVIGGMLVAEAFEHTPGIQTASLGTYLKTNLFLFLFAVGFYLLLRVLNIDLLWSVPIAKKWCANPDWIHIDTTPFAGLVRNLGVLFGLGFAINSEMFLLSCRGGNNYTLSFRLLCALTSLTILQLYHFLQIPTHEEHLFYVLSFCKSASIPLTVVAFIPYSVHMLMKQSGKKSQ


Regarding bad hits proteins, I’ve never personally looked specifically into these. But knowing that this particular MHC molecule prefers much hydrophobic peptides, I’ve checked on some hydrophilic bacterial proteins with SYFPEITHI (the current prediction algorithm that I’m using) and am sure that they give bad hits. So here they are.

Bad hits
> Leichmania Major hydrophilic surface protein 1
ETNQGGNASGSRKSAGGRTNEYDPKDDGFTPNNEDHCPKEDGHTGKNDDGGPKEDGHAPKNDDHAPKEDGHAAKNDDHAPKEDGHAQKNDGDVQ

> Leichmania Mexicana Hydrophilic acylated surface protein B
MGSSCTKDSTKEPEKRADNIDTTTGSNKKDGGHDHHQRTDGDGEKNDHDGEKADGDAGKHDDDHHQKTDGDGEKNDHDGEKADGDAGKHDDDHHQKTDGDGEKNDHDGEKADGDAGKHDDDDHHQKTDGDGEKNDHDGEKADGDAGNNEGEHHVGEGGDGDGAGNGNHPKRDTAAN


I hope this works well for you. Let me know should you need more info.

Thanks
Kailin


----



Hi Bosco,

The number of peptides for my dataset are as follows:
8mers - 21 peptides
9mers - 82 peptides
10mers - 36 peptides

As you can see the numbers are much lower than Kailin's so it will be interesting to see how this affects the output.

I've also attached three other files.
1) The bat proteome amino acid frequencies
2) The Hendra virus proteome
3) The bat proteome

The first should be used as background frequencies. The last two are proteomes I would like to search through for peptides. In particular, file 2, but if you can do file 3 that would be great too.

Hope that makes sense.

Cheers,
Amanda.