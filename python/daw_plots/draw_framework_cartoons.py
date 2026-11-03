"""Draw editable scientific SVG cartoons; optionally render using Inkscape.

Run from any directory: python3 python/draw_framework_cartoons.py --render
Only uses the Python standard library. No model calculations are performed.
"""
from pathlib import Path
from html import escape
import argparse
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'Results/ResultsFramework/cartoons'
INK = '#253538'
MUTED = '#657477'
POL = '#96b8b4'
POL_EDGE = '#456b68'
RNA = '#c58a25'
PAS = '#842f55'
BLUE = '#4b8ca9'
BLUE_EDGE = '#285e76'
GREEN = '#5a9b72'
GREY = '#737c80'


class Drawing:
    def __init__(self, title, height=710):
        self.parts = [f'''<svg xmlns="http://www.w3.org/2000/svg" width="180mm" height="{180*height/1560:.2f}mm" viewBox="0 0 1560 {height}">
<title>{escape(title)}</title>
<desc>Conceptual scientific cartoon. Molecule counts and spatial geometry are illustrative, not quantitative data.</desc>
<defs><marker id="arrow" viewBox="0 0 10 10" refX="8.5" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M0 0 10 5 0 10 2.5 5Z" fill="{MUTED}"/></marker></defs>
<style>text{{font-family:Arial,Helvetica,sans-serif;fill:{INK}}} path,line,polyline,circle,ellipse{{stroke-linecap:round;stroke-linejoin:round}}</style>
<rect width="1560" height="{height}" fill="white"/>''']
        self.text(30,48,'A',26,bold=True)
        self.text(74,48,title,26,bold=True)

    def add(self, svg):
        self.parts.append(svg)

    def text(self,x,y,value,size=22,color=INK,anchor='start',bold=False):
        self.add(f'<text x="{x}" y="{y}" font-size="{size}" fill="{color}" style="fill:{color}" text-anchor="{anchor}" font-weight="{600 if bold else 400}">{escape(value)}</text>')

    def path(self,d,color=INK,width=2,fill='none',extra=''):
        self.add(f'<path d="{d}" stroke="{color}" stroke-width="{width}" fill="{fill}" {extra}/>')

    def line(self,x1,y1,x2,y2,color=INK,width=2,extra=''):
        self.path(f'M{x1} {y1} L{x2} {y2}',color,width,extra=extra)

    def arrow(self,x1,y1,x2,y2,both=False):
        self.line(x1,y1,x2,y2,MUTED,2,'marker-end="url(#arrow)"'+(' marker-start="url(#arrow)"' if both else ''))

    def circle(self,x,y,r,fill,stroke='none',width=1):
        self.add(f'<circle cx="{x}" cy="{y}" r="{r}" fill="{fill}" stroke="{stroke}" stroke-width="{width}"/>')

    def group(self,name,transform=''):
        self.add(f'<g id="{name}" transform="{transform}">')

    def end(self):
        self.add('</g>')

    def factor(self,x,y,scale=1):
        self.add(f'<g transform="translate({x} {y}) scale({scale})">')
        self.path('M-17 -6 C-18 -18 -5 -23 6 -18 C18 -20 25 -6 18 5 C18 18 5 23 -6 17 C-20 17 -24 3 -17 -6Z',BLUE_EDGE,1.7,BLUE)
        self.text(0,7,'E',19,'white','middle')
        self.end()

    def site(self,x,y,formed=True):
        # One abstract site; do not imply two Ser2P marks per site.
        self.circle(x,y,6,GREEN if formed else '#d9dddf','white',1)

    def rate(self,x,y,sub,size=22):
        self.add(f'<text x="{x}" y="{y}" font-size="{size}" text-anchor="middle">'
                 '<tspan font-style="italic">k</tspan>'
                 f'<tspan baseline-shift="sub" font-size="{size*.7}">{escape(sub)}</tspan></text>')

    def polymerase(self,x,y,formed=5,bound=(1,3,4),scale=1,rna=True,
                   recognized=False,committed=False,rna_length='medium'):
        if committed and not recognized:
            raise ValueError('Commitment is shown only for a PAS-recognized complex')
        self.add(f'<g transform="translate({x} {y}) scale({scale})">')
        if committed:
            self.add(f'<ellipse cx="-29" cy="-122" rx="47" ry="33" fill="#f5eaf0" stroke="{PAS}" stroke-width="1.8"/>')
        self.path('M-50 -42 C-56 -61 -52 -85 -56 -100 C-60 -120 -68 -137 -81 -145 C-96 -159 -110 -166 -126 -165 C-142 -166 -160 -165 -171 -158 C-188 -150 -196 -137 -206 -130 C-222 -119 -236 -126 -236 -138',GREY,4.5)
        coords=[(-206,-130),(-171,-158),(-126,-165),(-81,-145),(-56,-100)]
        for i,(sx,sy) in enumerate(coords):
            self.site(sx,sy,i>=5-formed)
            if i in bound:
                if i<5-formed:
                    raise ValueError('An E factor must occupy a formed site')
                self.factor(sx,sy-20)
        self.path('M-71 5 C-83 -13 -76 -38 -60 -50 C-49 -67 -27 -62 -15 -55 C3 -69 28 -59 37 -46 C63 -42 78 -19 72 4 C78 27 54 48 29 49 C8 61 -17 53 -31 47 C-54 49 -69 31 -71 5Z',POL_EDGE,2.3,POL)
        self.text(0,7,'Pol II',22,INK,'middle')
        if rna:
            # The growing 3-prime end meets the upper Pol II outline, at Bezier t=.5.
            if recognized:
                self.path('M14.375 -60.625 C23 -67 25 -73 22 -77 C18 -93 3 -100 -13 -103 C-34 -106 -49 -119 -66 -114 C-96 -106 -105 -78 -133 -77 C-161 -75 -175 -89 -196 -83',RNA,6)
                self.text(-216,-80,'5′',19,RNA)
            elif rna_length=='long':
                self.path('M14.375 -60.625 C24 -76 25 -97 13 -111 C-5 -132 -30 -117 -42 -101 C-62 -74 -89 -75 -116 -83 C-146 -92 -162 -74 -185 -79 C-207 -85 -218 -77 -228 -74',RNA,6)
                self.text(-251,-71,'5′',19,RNA)
            else:
                self.path('M14.375 -60.625 C5 -71 6 -79 14 -90 C23 -115 5 -133 -11 -123 C-29 -109 -36 -83 -65 -77 C-97 -69 -122 -87 -143 -78',RNA,6)
                self.text(-166,-74,'5′',19,RNA)
            self.text(35,-64,'3′',19,RNA)
        if recognized:
            self.path('M-13 -103 C-23 -105 -31 -107 -39 -111',PAS,7)
            self.path('M-32 -121 Q-23 -125 -14 -117 L-8 -105 -19 -94 -32 -101 -38 -111Z',PAS,1.7,'#eedde5')
            self.text(-24,-105,'H',19,PAS,'middle')
        if committed:
            self.add('<g transform="translate(-2 -130) scale(.65)">')
            self.path('M-18 -6 C-21 -17 -7 -23 4 -19 C17 -23 25 -9 19 3 C22 16 8 22 -3 18 C-16 22 -25 8 -18 -6Z','#71478a',1.7,'#a579b5')
            self.text(0,7,'L',20,'white','middle')
            self.end()
        self.end()

    def dna(self,start,end,y):
        self.line(start,y,end,y,INK,2.7)
        self.line(start,y+12,end,y+12,INK,2.7)

    def tick(self,x,y,label,color=INK):
        self.line(x,y-9,x,y+23,color,3)
        self.text(x,y+49,label,20,color,'middle')

    def footer(self,y,lines):
        self.line(45,y,1515,y,'#dce3e2',1.1)
        for i,value in enumerate(lines):
            self.text(45,y+37+i*29,value,18,MUTED)

    def save(self,name):
        target=OUT/f'{name}.svg'
        target.write_text('\n'.join(self.parts+['</svg>']))
        return target


def resource_cartoon():
    d=Drawing('Local E loading draws on a shared factor pool',600)
    d.group('local-loading')
    d.text(60,115,'Loading on one gene',25,bold=True)
    d.dna(70,630,423)
    d.tick(85,423,'TSS');d.tick(565,423,'PAS',PAS)
    d.polymerase(405,398,scale=1.12)
    d.factor(433,210,.9)
    d.text(433,175,'Free E',20,BLUE_EDGE,'middle')
    d.arrow(401,198,348,198)
    d.arrow(348,224,401,224)
    d.rate(374,182,'Eon',21)
    d.rate(397,246,'Eoff',21)
    d.text(95,240,'CTD',21,MUTED)
    d.line(139,244,150,260,MUTED,1.2)
    d.arrow(480,330,568,330)
    d.text(72,525,'DNA',20)
    d.end()

    d.group('shared-pool')
    d.text(755,115,'Exchange across the gene population',25,bold=True)
    # This outline denotes an accounting pool, not a membrane-bound compartment.
    d.path('M940 182 C1010 154 1240 150 1321 186 C1382 215 1368 280 1304 300 C1218 325 1010 323 945 295 C883 270 883 207 940 182Z','#c6dce2',1.5,'#f0f6f7')
    d.text(1134,201,'Shared free E',23,BLUE_EDGE,'middle')
    for x,y in [(974,244),(1044,268),(1115,241),(1188,268),(1265,241)]:
        d.factor(x,y,.8)
    d.arrow(954,321,862,391,True)
    d.arrow(1128,334,1128,392,True)
    d.arrow(1315,321,1400,391,True)
    for i,(x,bound) in enumerate([(916,(3,)),(1181,(1,4)),(1455,(4,))]):
        d.dna(x-173,x+48,500)
        d.polymerase(x,488,bound=bound,scale=.56,rna=False)
        d.text(x-60,549,f'Gene {i+1}',20,INK,'middle')
    d.end()
    return d.save('Fig4A_resource_cartoon')


def length_cartoon():
    d=Drawing('Gene length sets the opportunity for E preloading',925)
    d.text(780,102,'Same free E and free Pol II for both gene lengths',23,BLUE_EDGE,'middle')
    for name,y,start,formed,bound in [('Short gene',360,665,1,()),('Long gene',710,170,5,(1,3,4))]:
        d.group(name.lower().replace(' ','-'))
        d.text(60,y-120,name,25,bold=True)
        d.dna(start,1230,y)
        d.tick(start,y,'TSS');d.tick(1115,y,'PAS',PAS)
        d.polymerase(1000,y-22,formed,bound,rna_length='long' if name=='Long gene' else 'medium')
        d.text(1200,y+49,'DNA',20)
        # Distance is shown explicitly; both polymerases are just upstream of the PAS.
        d.arrow(start+8,y+78,1107,y+78,True)
        d.text((start+1115)/2,y+108,'TSS-to-PAS distance',19,MUTED,'middle')
        if name=='Short gene':
            d.text(1254,y-55,'Less E loading',22,bold=True)
        else:
            d.text(1254,y-55,'More E loading',22,bold=True)
        d.end()
    d.group('site-legend')
    d.line(45,854,1515,854,'#dce3e2',1.1)
    d.site(65,891);d.text(87,898,'Ser2P site',19,GREEN)
    d.site(424,891,False);d.text(446,898,'Unphosphorylated',19,MUTED)
    d.factor(788,891,.65);d.text(816,898,'Bound E',19)
    d.text(1040,898,'M = 5',19,MUTED)
    d.end()
    return d.save('Fig5A_gene_length_cartoon')


def apa_cartoon():
    d=Drawing('PAS spacing defines a proximal commitment window',645)
    d.group('commitment-window')
    d.add('<rect x="255" y="163" width="1115" height="326" rx="14" fill="#f6f0f3"/>')
    d.text(780,119,'Opportunity to commit before the distal PAS is transcribed',24,INK,'middle')
    d.arrow(258,153,1366,153)
    d.text(780,199,'Inter-PAS distance, d',21,PAS,'middle')
    d.line(255,217,255,465,'#c8b4bd',1.4,'stroke-dasharray="5 7"')
    d.line(1370,217,1370,465,'#c8b4bd',1.4,'stroke-dasharray="5 7"')
    d.dna(95,1465,450)
    d.tick(255,450,'Proximal PAS',PAS)
    d.tick(1370,450,'Distal PAS',PAS)
    d.text(95,499,'DNA',20)
    d.text(590,235,'PAS recognition',22,anchor='middle',bold=True)
    d.text(1000,235,'CPA commitment',22,anchor='middle',bold=True)
    d.polymerase(685,428,bound=(1,3,4),scale=.9,recognized=True)
    d.polymerase(1095,428,bound=(1,3,4),scale=.9,recognized=True,committed=True)
    d.arrow(783,400,957,400)
    d.rate(870,377,'c')
    d.text(1267,295,'Distal PAS',20,PAS,'middle')
    d.text(1267,324,'not yet transcribed',19,MUTED,'middle')
    d.end()

    d.group('approximation')
    d.add('<text x="780" y="568" font-size="27" text-anchor="middle">Estimated proximal usage ≈ <tspan font-style="italic">F</tspan><tspan baseline-shift="sub" font-size="19">commit</tspan><tspan>(d)</tspan></text>')
    d.text(780,606,'Cumulative commitment probability evaluated at the inter-PAS distance',21,MUTED,'middle')
    d.end()
    return d.save('Fig6A_APA_cartoon')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--render',action='store_true')
    parser.add_argument('--inkscape',default=shutil.which('inkscape') or '/Applications/Inkscape.app/Contents/MacOS/inkscape')
    args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    drawings=[resource_cartoon(),length_cartoon(),apa_cartoon()]
    if args.render:
        for svg in drawings:
            subprocess.run([args.inkscape,str(svg),'--export-type=png','--export-width=2340',f'--export-filename={svg.with_suffix(".png")}'],check=True)
    for drawing in drawings:
        print(drawing.relative_to(ROOT))


if __name__=='__main__':
    main()
