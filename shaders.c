//
// Created by bruno on 2/14/17.
//

#define STRINGIFY(A) #A

/* vertex shader */
const char * depth_vertex_shader = STRINGIFY(
    uniform float point_radius;
    uniform float point_scale;
    uniform float density_scale;
    uniform float density_offset;
    uniform float max_z;
    uniform float min_z;
    vec3 vertex_position;
    varying vec3 pos_eye;
    void main(){

        pos_eye = vec3(gl_ModelViewMatrix * vec4(gl_Vertex.xyz, 1.0));
        float dist = length(pos_eye);

        gl_PointSize = point_radius * (point_scale / dist);
        gl_TexCoord[0] = gl_MultiTexCoord0;
        gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_Vertex.xyz, 1.0);

        gl_FrontColor = gl_Color;
    }
);

const char * depth_fragment_shader = STRINGIFY(
    uniform float point_radius;
    uniform float near;
    uniform float far;
    varying vec3 pos_eye;
    uniform vec3 light_dir = vec3(0.577, 0.577, 0.577);
    void main(){

        vec3 N;
        N.xy = gl_TexCoord[0].xy * vec2(2.0, -2.0) + vec2(-1.0, 1.0);
        float mag = dot(N.xy,N.xy);
        if(mag > 1.0f) discard;
        N.z = sqrt(1-mag);

        // point on surface on sphere in eye space
        vec4 sphere_pos_eye = vec4(pos_eye+N*point_radius, 1.0);
        vec4 clip_space_pos = gl_ProjectionMatrix * sphere_pos_eye;
        float norm_depth = clip_space_pos.z / clip_space_pos.w;


        float far = gl_DepthRange.far;
        float near = gl_DepthRange.near;
        gl_FragDepth=((far-near)*norm_depth+near+far)*0.5f;

        float diffuse = max(0.0, dot(light_dir, N));
        gl_FragColor = gl_Color*diffuse;

    }
);
const char * plane_vertex_shader = STRINGIFY(
    \n#version 400 core\n
    in vec3 position;

    uniform mat4 mView;
    uniform mat4 projection;

    out vec3 coord;

    void main() {

        gl_Position = projection * mView * vec4(position, 1.0);
        coord = position;
    }
);
const char * plane_fragment_shader = STRINGIFY(
    \n#version 400 core\n
    in vec3 coord;
    out vec3 color;

    void main() {
        vec3 col1 = vec3(0.7);
        vec3 col2 = vec3(0.9);

        color = mix(col1, col2, 0.5 * (mod(floor(coord.x*16) + floor(coord.y*16) + floor(coord.z*16), 2)));
    }
);

const char * water_depth_vertex_shader = STRINGIFY(
     \n#version 400 core\n

     in vec3 vertexPos;

     uniform mat4 projection;
     uniform mat4 mView;
     uniform vec2 screenSize;
     uniform float pointRadius;
     uniform float pointScale;

     out vec3 pos;

     void main() {
          vec4 viewPos = mView * vec4(vertexPos.xyz, 1.0);
          gl_Position = projection * viewPos;
          pos = viewPos.xyz;

          gl_PointSize = pointScale * (pointRadius / gl_Position.w);
     }
);

const char * water_depth_fragment_shader = STRINGIFY(
    \n#version 400 core\n

    in vec3 pos;

    uniform mat4 mView;
    uniform mat4 projection;
    uniform float pointRadius;
    uniform float pointScale;
    void main() {

        //calculate normal
        vec3 normal;

        normal.xy = gl_PointCoord * 2.0 - 1.0;
        float r2 = dot(normal.xy, normal.xy);

        if (r2 > 1.0) {
            discard;
        }

        normal.z = sqrt(1.0 - r2);

        //calculate depth
        vec4 pixelPos = vec4(pos + normal * pointRadius, 1.0);
        vec4 clipSpacePos = projection * pixelPos;

        float far = gl_DepthRange.far;
        float near = gl_DepthRange.near;
        gl_FragDepth=((far-near)*(clipSpacePos.z / clipSpacePos.w)+near+far)*0.5f;
        //gl_FragDepth = (clipSpacePos.z / clipSpacePos.w) * 0.5f + 0.5f;
}
);


/* BLUR SHADER */
const char * blur_vertex_shader = STRINGIFY(
    \n#version 400 core\n
    in vec2 vertexPos;
    out vec2 coord;
    void main(){
        coord = 0.5f * vertexPos + 0.5f;
        gl_Position = vec4(vertexPos, 0.0f, 1.0f);
    }
);

const char * blur_fragment_shader = STRINGIFY(
    \n#version 400 core\n
    in vec2 coord;

    uniform sampler2D depthMap;
    uniform vec2 screenSize;
    uniform mat4 projection;
    uniform vec2 blurDir;
    uniform float filterRadius;
    uniform float blurScale;

    const float blurDepthFalloff = 65.0f;

    void main() {
        float depth = texture(depthMap, coord).x;

        if (depth <= 0.0f) {
            gl_FragDepth = 0;
            return;
        }

        if (depth >= 1.0f) {
            gl_FragDepth = depth;
            return;
        }

        float sum = 0.0f;
        float wsum = 0.0f;

        for (float x = -filterRadius; x <= filterRadius; x += 1.0f) {
            float s = texture(depthMap, coord + x * blurDir).x;

            //if (s >= 1.0f) continue;

            float r = x * blurScale;
            float w = exp(-r*r);

            float r2 = (s - depth) * blurDepthFalloff;
            float g = exp(-r2*r2);

            sum += s * w * g;
            wsum += w * g;
        }

        if (wsum > 0.0f) {
            sum /= wsum;
        }

        gl_FragDepth = sum;
    }
);

/* thickness frag */
const char * thickness_fragment_shader = STRINGIFY(

     \n#version 400 core\n
     in vec3 pos;
     uniform mat4 mView;
     uniform mat4 projection;

     out float thickness;

     void main() {
        //calculate normal
        vec3 normal;
        normal.xy = gl_PointCoord * 2.0 - 1.0;
        float r2 = dot(normal.xy, normal.xy);

        if (r2 > 1.0f) {
            discard;
        }

        normal.z = sqrt(1 - r2);

        thickness = normal.z * 0.07f;

    }
);

const char * shadow_depth_vertex_shader = STRINGIFY(
        \nversion 400 core\n

        void main(){

}
);

const char * shadow_depth_fragment_shader = STRINGIFY(
        \nversion 400 core\n

        void main(){

}
);

const char * shadow_map_vertex_shader = STRINGIFY(
    \nversion 400 core\n

    void main(){

    }
);

const char * shadow_map_fragment_shader = STRINGIFY(
        \nversion 400 core\n

        void main(){

}
);

/* fluid final shader */
const char * fluidfinal_fragment_shader = STRINGIFY(
        \n#version 400 core\n

        in vec2 coord;

        uniform vec4 color;
        uniform sampler2D depthMap;
        uniform sampler2D thicknessMap;
        uniform sampler2D sceneMap;
        uniform mat4 projection;
        uniform mat4 mView;
        uniform vec2 invTexScale;

        out vec4 fragColor;

        const vec3 lightDir = vec3(0, 1, 1);
        //const vec3 lightDir = vec3(0.577, 0.577, 0.577);
        const vec3 lightPos = vec3(0, 1000, 1000);
        const float shininess = 500.0;
        const float fresPower = 25.0f;
        const float fresScale = 0.9;
        const float fresBias = 0.1;

        vec3 uvToEye(vec2 p, float z) {
             vec2 pos = p * 2.0f - 1.0f;
             vec4 clipPos = vec4(pos, z, 1.0f);
             vec4 viewPos = inverse(projection) * clipPos;
             return viewPos.xyz / viewPos.w;
        }

        void main() {
           vec4 scene = texture(sceneMap, coord);
           float depth = texture(depthMap, coord).x;

           if (depth == 0.0f) {
               fragColor = vec4(0);
               return;
           }

           if (depth == 1.0) {
              fragColor = scene;
              return;
           }

           // reconstruct eye space pos from depth
           vec3 eyePos = uvToEye(coord, depth);

           // finite difference approx for normals, can't take dFdx because
           // the one-sided difference is incorrect at shape boundaries
           vec3 zl = eyePos - uvToEye(coord - vec2(invTexScale.x, 0.0), texture(depthMap, coord - vec2(invTexScale.x, 0.0)).x);
           vec3 zr = uvToEye(coord + vec2(invTexScale.x, 0.0), texture(depthMap, coord + vec2(invTexScale.x, 0.0)).x) - eyePos;
           vec3 zt = uvToEye(coord + vec2(0.0, invTexScale.y), texture(depthMap, coord + vec2(0.0, invTexScale.y)).x) - eyePos;
           vec3 zb = eyePos - uvToEye(coord - vec2(0.0, invTexScale.y), texture(depthMap, coord - vec2(0.0, invTexScale.y)).x);

           vec3 dx = zl;
           vec3 dy = zt;

           if (abs(zr.z) < abs(zl.z))
               dx = zr;

           if (abs(zb.z) < abs(zt.z))
               dy = zb;

           vec3 normal = normalize(cross(dx, dy));
           vec4 worldPos = inverse(mView) * vec4(eyePos, 1.0);

           //Phong specular
           vec3 l = (mView * vec4(lightDir, 0.0)).xyz;
           vec3 viewDir = -normalize(eyePos);
           vec3 halfVec = normalize(viewDir + l);
           float specular = pow(max(0.0f, dot(normal, halfVec)), shininess);

           vec2 texScale = vec2(0.75, 1.0);
           float refractScale = 1.33f * 0.025f;
           refractScale *= smoothstep(0.1, 0.4, worldPos.y);
           vec2 refractCoord = coord + normal.xy*refractScale*texScale;

           //float thickness = max(texture(thicknessMap, refractCoord).x, 0.3);
           float thickness = max(texture(thicknessMap, coord).x, 0.3);
           vec3 transmission = exp(-(vec3(1.0)-color.xyz)*thickness);
           //vec3 transmission = (1.0-(1.0-color.xyz)*thickness*0.8)*color.w;

           vec3 refract = texture(sceneMap, refractCoord).xyz*transmission;

           vec3 lVec = normalize(worldPos.xyz-lightPos);
           float attenuation = max(smoothstep(0.95, 1.0, abs(dot(lVec, -lightDir))), 0.05);
           float ln = dot(l, normal)*attenuation;

           //Fresnel
           float fresnel = fresBias + fresScale * pow(1.0f - max(dot(normal, viewDir), 0.0), fresPower);

           //Diffuse light
           vec3 diffuse = color.xyz * mix(vec3(0.29, 0.379, 0.59), vec3(1.0), (ln*0.5 + 0.5)) * (1 - color.w);
           //vec3 diffuse = color.xyz * mix(vec3(0.29, 0.379, 0.59), vec3(1.0), (ln*0.5 + 0.5));

           vec3 skyColor = vec3(0.1, 0.2, 0.4)*1.2;
           vec3 groundColor = vec3(0.1, 0.1, 0.2);
           vec3 rEye = reflect(viewDir, normal).xyz;
           vec3 rWorld = (inverse(mView)*vec4(rEye, 0.0)).xyz;
           vec3 reflect = vec3(1.0) + mix(groundColor, skyColor, smoothstep(0.15, 0.25, rWorld.y));

           //Compositing everything
           vec3 finalColor = diffuse + (mix(refract, reflect, fresnel) + specular) * color.w;

           fragColor = vec4(finalColor, 1.0);
           gl_FragDepth = depth;
    }
);

const char * fluidfinal_vertex_shader = STRINGIFY(
    \n#version 400 core\n
    in vec2 vertexPos;
    out vec2 coord;

    void main() {
        coord = 0.5f * vertexPos + 0.5f;
        gl_Position = vec4(vertexPos, 0.0f, 1.0);
    }
);


const char * final_fragment_shader = STRINGIFY(
    \n#version 400 core\n

    in vec2 coord;
    uniform sampler2D fluidMap;
    out vec4 fragColor;

    void main() {
      
         vec4 fluid = texture(fluidMap, coord).xyzw;
         fragColor = fluid;
    }
);

const char * final_vertex_shader = STRINGIFY(
    \n#version 400 core\n
    in vec2 vertexPos;

    out vec2 coord;

    void main() {
         coord = 0.5f * vertexPos + 0.5f;
         gl_Position = vec4(vertexPos, 0.0f, 1.0);
    }
);
