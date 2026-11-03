"""Editable review draft of Figure 2A. Does not modify manuscript documents.

User clarification: kc is commitment; subsequent cleavage is not modeled.
Run: python3 python/draw_figure2a_revision.py
Render the resulting SVG with Inkscape for a PNG preview.
"""
from pathlib import Path
from python.daw_plots.draw_framework_cartoons import (
    Drawing, INK, MUTED, POL, POL_EDGE, RNA, PAS, BLUE, BLUE_EDGE, GREEN, GREY,
)

OUT = Path(__file__).resolve().parents[1] / 'Results/ResultsFramework/figure2a_revision'


class Revision(Drawing):
    def site(self, x, y, formed=True):
        # One abstract site symbol; no definition in terms of paired Ser2P marks.
        self.circle(x, y, 6, GREEN if formed else '#d9dddf', 'white', 1)

    def rate(self, x, y, sub, prefix='', size=22, anchor='middle'):
        self.add(f'<text x="{x}" y="{y}" font-size="{size}" text-anchor="{anchor}">'
                 f'<tspan font-style="italic">{prefix}k</tspan>'
                 f'<tspan baseline-shift="sub" font-size="{size*.7}">{sub}</tspan></text>')

    def polymerase(self, x, y, formed=5, bound=(1,3,4), scale=1,
                   rna=True, recognized=False, rna_length='medium'):
        self.group(f'pol-{len(self.parts)}', f'translate({x} {y}) scale({scale})')
        super().polymerase(0,0,formed=formed,bound=bound,rna=False,recognized=False)
        if rna:
            # (14.375,-60.625) lies on the upper Pol II outline (Bezier t=.5).
            # The visible RNA ends exactly at that edge rather than inside the body.
            if recognized:
                curve='M14.375 -60.625 C23 -67 25 -73 22 -77 C18 -93 3 -100 -13 -103 C-34 -106 -49 -119 -66 -114 C-96 -106 -105 -78 -133 -77 C-161 -75 -175 -89 -196 -83'
                five=(-216,-80)
            elif rna_length == 'short':
                curve='M14.375 -60.625 C8 -77 -8 -83 -27 -81'
                five=(-51,-79)
            elif rna_length == 'long':
                curve='M14.375 -60.625 C24 -76 25 -97 13 -111 C-5 -132 -30 -117 -42 -101 C-62 -74 -89 -75 -116 -83 C-146 -92 -162 -74 -185 -79 C-207 -85 -218 -77 -228 -74'
                five=(-251,-71)
            else:
                curve='M14.375 -60.625 C5 -71 6 -79 14 -90 C23 -115 5 -133 -11 -123 C-29 -109 -36 -83 -65 -77 C-97 -69 -122 -87 -143 -78'
                five=(-166,-74)
            self.path(curve,RNA,6)
            self.text(*five,'5′',19,RNA)
            self.text(35,-64,'3′',19,RNA)
        if recognized:
            self.path('M-13 -103 C-23 -105 -31 -107 -39 -111',PAS,7)
            self.path('M-32 -121 Q-23 -125 -14 -117 L-8 -105 -19 -94 -32 -101 -38 -111Z',PAS,1.7,'#eedde5')
            self.text(-24,-105,'H',19,PAS,'middle')
        self.end()

    def molecule(self, x, recognized=False, committed=False):
        self.group(f'molecule-{x}', f'translate({x} 365) scale(.88)')
        if committed:
            self.add(f'<ellipse cx="-29" cy="-122" rx="47" ry="33" fill="#f5eaf0" stroke="{PAS}" stroke-width="1.8"/>')
        self.polymerase(0, 0, formed=5, bound=(1,3,4), recognized=recognized)
        if committed:
            self.group('commitment-L','translate(-2 -130) scale(.65)')
            self.path('M-18 -6 C-21 -17 -7 -23 4 -19 C17 -23 25 -9 19 3 C22 16 8 22 -3 18 C-16 22 -25 8 -18 -6Z','#71478a',1.7,'#a579b5')
            self.text(0,7,'L',20,'white','middle')
            self.end()
        self.end()


def draw():
    d = Revision('CTD loading and PAS recognition lead to CPA commitment', 890)
    # Main sequence. All molecular objects and labels remain editable SVG.
    for x, title in [(45,'1  E loading'), (425,'2  PAS recognition'),
                     (805,'3  CPA commitment'), (1185,'4  RNA cleavage')]:
        d.text(x,113,title,25,bold=True)
    d.text(45,145,'During elongation',20,MUTED)
    d.text(425,145,'Reversible',20,MUTED)
    d.text(805,145,'Irreversible',20,PAS)
    d.text(1185,145,'Downstream biological event',20,MUTED)

    for start in [50,430,810,1190]:
        d.dna(start,start+302,389)
    d.molecule(276)
    d.molecule(656,recognized=True)
    d.molecule(1036,recognized=True,committed=True)

    # E exchange is shown in the expanded mid-gene snapshot below.
    d.text(48,202,'CTD',21,MUTED)
    d.line(86,208,98,235,MUTED,1.2)
    d.text(50,434,'DNA',20)
    d.text(145,324,'RNA',20,RNA)
    d.text(201,471,'E-loaded Pol II',22,anchor='middle')

    d.text(715,288,'PAS',20,PAS)
    d.line(707,282,653,276,MUTED,1.2)
    d.text(581,471,'Recognition complex',22,anchor='middle')
    d.text(961,471,'Committed to CPA',22,anchor='middle')

    # Recognition is reversible; kc ends the modeled sequence at commitment.
    d.arrow(360,343,418,343)
    d.arrow(418,361,360,361)
    d.arrow(740,351,798,351)
    d.rate(769,326,'c')

    # Dashed continuation distinguishes the unmodeled downstream illustration.
    d.line(1120,351,1178,351,MUTED,2,'stroke-dasharray="6 6" marker-end="url(#arrow)"')
    d.line(1149,165,1149,487,'#d3dada',1.3,'stroke-dasharray="4 6"')

    # Cleaved upstream RNA with a PAS and a visibly separate downstream fragment.
    d.path('M1218 257 C1244 245 1265 269 1294 257 C1311 249 1318 233 1341 241 C1353 246 1365 251 1375 251',RNA,5.3)
    d.path('M1341 241 Q1355 248 1366 250',PAS,6.3)
    d.text(1192,261,'5′',19,RNA)
    d.text(1289,208,'Cleaved RNA',21,RNA,'middle')
    d.circle(1390,218,6,'white',INK,1.7)
    d.circle(1409,218,6,'white',INK,1.7)
    d.path('M1394 223 L1411 256 M1405 223 L1389 256',INK,1.8)
    # Pol II remains as biological context: no claim that cleavage immediately terminates it.
    d.group('post-cleavage-polymerase','translate(1420 365) scale(.88)')
    # Omit CTD detail here rather than inventing post-commitment site/E states.
    d.path('M-71 5 C-83 -13 -76 -38 -60 -50 C-49 -67 -27 -62 -15 -55 C3 -69 28 -59 37 -46 C63 -42 78 -19 72 4 C78 27 54 48 29 49 C8 61 -17 53 -31 47 C-54 49 -69 31 -71 5Z',POL_EDGE,2.3,POL)
    d.text(0,7,'Pol II',22,INK,'middle')
    d.end()
    d.path('M1420 255 C1451 261 1456 286 1432.65 311.65',RNA,5.3)
    d.text(1451,309,'3′',17,RNA)

    # Expanded step 1: three illustrative snapshots on the same gene.
    d.group('loading-detail-panel','translate(0 -65)')
    d.line(45,571,1515,571,'#dce3e2',1.1)
    d.text(45,610,'Step 1 in detail: Ser2P phosphorylation and E loading during elongation',23,bold=True)
    d.group('loading-along-gene')
    d.dna(70,1515,879)
    snapshots=[(340,'At the TSS',0,(),'short'),
               (840,'Mid-gene',3,(3,),'medium'),
               (1340,'Near the PAS',5,(1,3,4),'long')]
    for x,label,formed,bound,length in snapshots:
        d.text(x-(110 if x==840 else 75),660,label,23,anchor='middle',bold=True)
        d.polymerase(x,855,formed=formed,bound=bound,scale=.9,rna_length=length)
    # Show one phosphorylated CTD segment as the alternative state of this site.
    d.group('tss-phosphorylation-cycle')
    d.path('M381 718 Q400 708 419 718',GREY,3)
    d.site(400,713,True)
    d.path('M276 718 Q314 681 375 704',MUTED,1.7,extra='marker-end="url(#arrow)"')
    d.path('M378 734 Q324 759 276 733',MUTED,1.7,extra='marker-end="url(#arrow)"')
    d.rate(330,681,'Pon',size=21)
    d.rate(335,767,'Poff',size=21)
    d.end()
    d.factor(874,700,.9)
    d.text(874,665,'Free E',20,BLUE_EDGE,'middle')
    # Short, parallel exchange arrows have distinct vertical lanes.
    d.arrow(848,689,795,689)
    d.arrow(795,716,848,716)
    d.rate(818,674,'Eon',size=21)
    d.rate(818,733,'Eoff',size=21)
    d.arrow(453,788,567,788)
    d.arrow(953,788,1067,788)
    d.rate(510,773,'p',size=24)
    d.rate(1010,773,'p',size=24)
    # Place the TSS guide below the silhouette so it does not cut through Pol II.
    d.line(340,907,340,916,INK,2)
    d.text(340,929,'TSS',20,anchor='middle')
    d.text(70,929,'DNA',20)
    d.tick(1490,879,'PAS',PAS)
    d.circle(556,922,6,'#d9dddf','white',1)
    d.text(572,929,'Unphosphorylated',19,MUTED)
    d.circle(802,922,6,GREEN,'white',1)
    d.text(818,929,'Ser2P site',19,GREEN)
    d.factor(1011,921,.6)
    d.text(1033,929,'Bound E',19,BLUE_EDGE)
    d.end()
    d.end()

    OUT.mkdir(parents=True,exist_ok=True)
    output = OUT/'Fig2A_revision.svg'
    output.write_text('\n'.join(d.parts+['</svg>']))
    print(output)


if __name__ == '__main__':
    draw()
