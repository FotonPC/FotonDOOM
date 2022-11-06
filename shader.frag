
uniform vec2 u_resolution;
uniform vec2 u_mouse;
uniform vec3 u_pos;
uniform float u_time;
uniform sampler2D u_sample;
uniform sampler2D u_h_map;
uniform sampler2D u_flat_sprites_info;
uniform sampler2D u_textures[16];
uniform int rx;
uniform int ry;
uniform float u_viewnetka;
uniform float u_blur;
const vec4 bitEnc = vec4(1.,255.,65025.,16581375.);
const vec4 bitDec = 1./bitEnc;
vec4 EncodeFloatRGBA (float v) {
    vec4 enc = bitEnc * v;
    enc = fract(enc);
    enc -= enc.yzww * vec2(1./255., 0.).xxxy;
    return enc;
}
float DecodeFloatRGBA (vec4 v) {
    return dot(v, bitDec);
}

const float MAX_DIST = 99999.0;
const int MAX_REF = 8;
vec3 light = normalize(vec3(-0.5, 0.75, -1.0));
void main() {
	vec2 xy = gl_TexCoord[0].xy;
    vec4 texColor = texture(u_sample, xy);
    vec2 xy1 = xy;
    xy1 = xy1 - 1/u_resolution.xy*2;
    vec2 xy2 = xy;
    xy2 = xy2 + 1/u_resolution.xy*2;
    vec2 xy3 = xy;
    xy3.x = xy3.x - 1/u_resolution.x*2;
    xy3.y = xy3.y + 1/u_resolution.y*2;
    vec2 xy4 = xy;
    xy4.x = xy4.x + 1/u_resolution.x*2;
    xy4.y = xy4.y - 1/u_resolution.y*2;
    vec4 texColor1 = texture(u_sample, xy1);
    vec4 texColor2 = texture(u_sample, xy2);
    vec4 texColor3 = texture(u_sample, xy3);
    vec4 texColor4 = texture(u_sample, xy4);
    texColor = (texColor + texColor1 * u_blur * abs(xy.x - 0.5) * abs(xy.x - 0.5) +texColor2 * u_blur * abs(xy.x - 0.5) * abs(xy.x - 0.5)+texColor3 * u_blur * abs(xy.x - 0.5) * abs(xy.x - 0.5)+texColor4 * u_blur * abs(xy.x - 0.5) * abs(xy.x - 0.5))/(4* u_blur * abs(xy.x - 0.5) * abs(xy.x - 0.5) + 1);
    texColor.rgb = texColor.rgb * min(vec3(1), (0.5-(sqrt(abs(xy.x - 0.5)*abs(xy.x - 0.5) + abs(xy.y-0.5)*abs(xy.y-0.5))) + u_viewnetka));
    vec4 h = texture(u_h_map, xy);
    texColor.rgb = texColor.rgb*(1-h.r);
    gl_FragColor = texColor;
} ;