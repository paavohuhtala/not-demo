
in vec2 v_uv;

uniform float time;

out vec4 frag;

void main() {
    frag = vec4(v_uv.x, 0.0, v_uv.y * abs(sin(time)), 1.0);
}
