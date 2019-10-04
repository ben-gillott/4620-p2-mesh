import sys

head_fmt = """\
<?xml version="1.0" encoding="UTF-8" ?>
<!-- 
{title:s}
(This file is generated by make_scenes.py. Edit and re-run
that script rather than editing this file.)
  -->
<scene>
"""

tail_fmt = """\
</scene>
"""

def write_file(fname, title, parts):
    with open(fname, 'w') as f:
        f.write(head_fmt.format(title=title))
        for format, args in parts:
            f.write(format.format(**args))
        f.write(tail_fmt.format(title=title))

bunny_cam = """
  <exposure>20</exposure>
  <camera type="PerspectiveCamera">
    <viewPoint>4 6 8</viewPoint>
    <viewDir>-4 -6 -8</viewDir>
    <viewUp>0 1 0</viewUp>
    <projDistance>2</projDistance>
    <viewWidth>0.5</viewWidth>
    <viewHeight>0.5</viewHeight>
  </camera>
  <image>
    400 400
  </image>
"""

bunny_shaders_ufacet = """
  <shader name="ground" type="Lambertian">
    <diffuseColor>.04 0.12 0.024</diffuseColor>
  </shader>
  <shader name="bunny" type="Microfacet{variant:s}">
    <diffuseColor>0.8 0.4 0.1</diffuseColor>
    <nt>1.6</nt>
    <alpha>{alpha:f}</alpha>
  </shader>
"""

bunny_shaders_phong = """
  <shader name="ground" type="Lambertian">
    <diffuseColor>.04 0.12 0.024</diffuseColor>
  </shader>
  <shader name="bunny" type="Phong">
    <diffuseColor>0.8 0.4 0.1</diffuseColor>
    <specularColor>0.2 0.1 0.025</specularColor>
    <exponent>25</exponent>
  </shader>
"""

bunny_scene = """
  <surface type="Mesh">
    <shader ref="bunny" />
    <data>../../meshes/{mesh:s}.obj</data>
  </surface>
  <surface type="Box">
    <minpt>-3 -2 -3</minpt>
    <maxpt>3 -0.9 3</maxpt>
    <shader ref="ground" />
  </surface>
  
  <light>
    <position>3 10 5</position>
    <intensity>25 25 25</intensity>
  </light>
"""

write_file('bunny.xml', 'Faceted Stanford bunny', [
    (bunny_cam, {}),
    (bunny_shaders_ufacet, {
        'variant': 'Beckmann',
        'alpha': 0.2
    }),
    (bunny_scene, {
        'mesh': 'bunny'
    })
])

write_file('bunny_norms_beck1.xml', 
    'Smooth Stanford bunny with rougher Beckmann microfacet material', [
    (bunny_cam, {}),
    (bunny_shaders_ufacet, {
        'variant': 'Beckmann',
        'alpha': 0.15
    }),
    (bunny_scene, {
        'mesh': 'bunnyNV'
    })
])

write_file('bunny_norms_beck2.xml', 
    'Smooth Stanford bunny with smoother Beckmann microfacet material', [
    (bunny_cam, {}),
    (bunny_shaders_ufacet, {
        'variant': 'Beckmann',
        'alpha': 0.05
    }),
    (bunny_scene, {
        'mesh': 'bunnyNV'
    })
])

write_file('bunny_norms_ggx1.xml', 
    'Smooth Stanford bunny with rougher GGX microfacet material', [
    (bunny_cam, {}),
    (bunny_shaders_ufacet, {
        'variant': 'GGX',
        'alpha': 0.15
    }),
    (bunny_scene, {
        'mesh': 'bunnyNV'
    })
])

write_file('bunny_norms_ggx2.xml', 
    'Smooth Stanford bunny with smoother GGX microfacet material', [
    (bunny_cam, {}),
    (bunny_shaders_ufacet, {
        'variant': 'GGX',
        'alpha': 0.05
    }),
    (bunny_scene, {
        'mesh': 'bunnyNV'
    })
])

write_file('bunny_Phong.xml', 
    'Smooth Stanford bunny with Phong material', [
    (bunny_cam, {}),
    (bunny_shaders_phong, {}),
    (bunny_scene, {
        'mesh': 'bunny'
    })
])

write_file('bunny_norms_Phong.xml', 
    'Smooth Stanford bunny with Phong material', [
    (bunny_cam, {}),
    (bunny_shaders_phong, {}),
    (bunny_scene, {
        'mesh': 'bunnyNV'
    })
])

teapot_cam = """
  <exposure>15</exposure>
  <camera type="PerspectiveCamera">
    <viewPoint>4 6 6</viewPoint>
    <viewDir>-3.5 -4.65 -6</viewDir>
    <viewUp>0 1 0</viewUp>
    <projDistance>5</projDistance>
    <viewWidth>5</viewWidth>
    <viewHeight>3</viewHeight>
  </camera>
  <image>
    400 240
  </image>
"""

teapot_shaders = """
  <shader name="teapot" type="Microfacet{variant:s}">
    <diffuseColor>.1 .15 .4</diffuseColor>
    <alpha>{alpha:f}</alpha>
    <nt>1.5</nt>
    <mirrorCoefficient>0.3 0.3 0.3</mirrorCoefficient>
  </shader>
"""

teapot_scene = """
  <surface type="Mesh">
    <shader ref="teapot" />
    <data>../../meshes/teapot.obj</data>
  </surface>
  <light>
    <position>7 6 10</position>
    <intensity>40 40 40</intensity>
  </light>
  <light>
    <position>4 3 -6</position>
    <intensity>30 30 30</intensity>
  </light>
"""

write_file('teapot.xml', 
    'Utah teapot with smoother Beckmann microfacet material', [
    (teapot_cam, {}),
    (teapot_shaders, {
        'variant': 'Beckmann',
        'alpha': 0.1
    }),
    (teapot_scene, {})
])

write_file('teapot2.xml', 
    'Utah teapot with rougher Beckmann microfacet material', [
    (teapot_cam, {}),
    (teapot_shaders, {
        'variant': 'Beckmann',
        'alpha': 0.2
    }),
    (teapot_scene, {})
])

write_file('teapot-GGX.xml', 
    'Utah teapot with smoother GGX microfacet material', [
    (teapot_cam, {}),
    (teapot_shaders, {
        'variant': 'GGX',
        'alpha': 0.1
    }),
    (teapot_scene, {})
])

write_file('teapot2-GGX.xml', 
    'Utah teapot with rougher GGX microfacet material', [
    (teapot_cam, {}),
    (teapot_shaders, {
        'variant': 'GGX',
        'alpha': 0.2
    }),
    (teapot_scene, {})
])

