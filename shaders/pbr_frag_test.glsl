varying vec3 vViewPosition;
varying vec3 vNormal;

varying vec2 vUv;

// uniforms
uniform float metallic;
uniform vec3 albedo;
uniform int diffuseType;

uniform sampler2D roughnessMap;
uniform sampler2D diffuseMap;

// defines
#define PI 3.14159265359
#define PI2 6.28318530718
#define RECIPROCAL_PI 0.31830988618
#define RECIPROCAL_PI2 0.15915494
#define LOG2 1.442695
#define EPSILON 1e-6

struct IncidentLight {
  vec3 color;
  vec3 direction;
  bool visible;
};

struct ReflectedLight {
  vec3 directDiffuse;
  vec3 directSpecular;
  vec3 indirectDiffuse;
  vec3 indirectSpecular;
};

struct GeometricContext {
  vec3 position;
  vec3 normal;
  vec3 viewDir;
};

struct Material {
  vec3 diffuseColor;
  float specularRoughness;
  vec3 specularColor;
};

// lights

bool testLightInRange(const in float lightDistance, const in float cutoffDistance) {
  return any(bvec2(cutoffDistance == 0.0, lightDistance < cutoffDistance));
}

float punctualLightIntensityToIrradianceFactor(const in float lightDistance, const in float cutoffDistance, const in float decayExponent) {
  if (decayExponent > 0.0) {
    return pow(saturate(-lightDistance / cutoffDistance + 1.0), decayExponent);
  }
  
  return 1.0;
}

struct DirectionalLight {
  vec3 direction;
  vec3 color;
};

void getDirectionalDirectLightIrradiance(const in DirectionalLight directionalLight, const in GeometricContext geometry, out IncidentLight directLight) {
  directLight.color = directionalLight.color;
  directLight.direction = directionalLight.direction;
  directLight.visible = true;
}

struct PointLight {
  vec3 position;
  vec3 color;
  float distance;
  float decay;
};

void getPointDirectLightIrradiance(const in PointLight pointLight, const in GeometricContext geometry, out IncidentLight directLight) {
  vec3 L = pointLight.position - geometry.position;
  directLight.direction = normalize(L);
  
  float lightDistance = length(L);
  if (testLightInRange(lightDistance, pointLight.distance)) {
    directLight.color = pointLight.color;
    directLight.color *= punctualLightIntensityToIrradianceFactor(lightDistance, pointLight.distance, pointLight.decay);
    directLight.visible = true;
  } else {
    directLight.color = vec3(0.0);
    directLight.visible = false;
  }
}

struct SpotLight {
  vec3 position;
  vec3 direction;
  vec3 color;
  float distance;
  float decay;
  float coneCos;
  float penumbraCos;
};

void getSpotDirectLightIrradiance(const in SpotLight spotLight, const in GeometricContext geometry, out IncidentLight directLight) {
  vec3 L = spotLight.position - geometry.position;
  directLight.direction = normalize(L);
  
  float lightDistance = length(L);
  float angleCos = dot(directLight.direction, spotLight.direction);
  
  if (all(bvec2(angleCos > spotLight.coneCos, testLightInRange(lightDistance, spotLight.distance)))) {
    float spotEffect = smoothstep(spotLight.coneCos, spotLight.penumbraCos, angleCos);
    directLight.color = spotLight.color;
    directLight.color *= spotEffect * punctualLightIntensityToIrradianceFactor(lightDistance, spotLight.distance, spotLight.decay);
    directLight.visible = true;
  } else {
    directLight.color = vec3(0.0);
    directLight.visible = false;
  }
}

// light uniforms
#define LIGHT_MAX 4
uniform DirectionalLight directionalLights[LIGHT_MAX];
uniform PointLight pointLights[LIGHT_MAX];
uniform SpotLight spotLights[LIGHT_MAX];
uniform int numDirectionalLights;
uniform int numPointLights;
uniform int numSpotLights;

// BRDFs

vec3 GGXOrenNayarDiffuse(vec3 color, float dotNL, float dotNV, float dotNH, float dotLH, float dotLV, float roughness, float a) {
  float f0 = 0.04;
  float theta_i = acos(dotNL);
  float theta_r = acos(dotNV);
  float cos_phi_diff = (dotLV - dotNL*dotNV) / (sin(theta_i)*sin(theta_r) + EPSILON);
  
  float Fr1 = (1.0 - (0.542026*a + 0.303573*roughness) / (a + 1.36053));
  float Fr2 = (1.0 - (pow(1.0 - dotNV, 5.0 - 4.0*a)) / (a + 1.36053));
  float Fr3 = (-0.733996*a*roughness + 1.50912*a - 1.16402*roughness);
  float Fr4 = (pow(1.0 - dotNV, 1.0 + (1.0 / (39.0*a*a + 1.0))));
  float Fr = Fr1*Fr2*(Fr3*Fr4+1.0);
  float Lm1 = (max(1.0 - (2.0*roughness), 0.0)*(1.0 - pow(1.0 - dotNL, 5.0)) + min(2.0*roughness, 1.0)); 
  float Lm2 = ((1.0 - 0.5*roughness)*dotNL + (0.5*roughness)*(pow(dotNL, 2.0)));
	float Lm = Lm1 * Lm2;
	float Vd1 = (a / ((a + 0.09)*(1.31072 + 0.995584*dotNV)));
  float Vd2 = (1.0 - (pow(1.0 - dotNL, (1.0 - 0.3726732*(dotNV*dotNV)) / (0.188566 + 0.38841*dotNV))));
  float Vd = Vd1*Vd2;
  float Bp = dotLV - (dotNV*dotNL);
	if (Bp < 0.0) Bp *= 1.4*dotNV*dotNL;
  float L1 = 1.05*(1.0 - f0)*(Fr*Lm + Vd*Bp);
  return color*(L1*RECIPROCAL_PI);
}

vec3 GGXApproxDiffuse(vec3 color, float dotNL, float dotNV, float dotNH, float dotLH, float dotLV, float roughness, float a) {
  float facing = 0.5 + 0.5*dotLV;
  float rough  = facing*(0.9-0.4*facing)*((0.5+dotNH)/dotNH);
  float Smooth = 1.05*(1.0-pow(1.0-dotNL, 5.0)) * (1.0-pow(1.0-dotNV, 5.0));
  float single = mix(Smooth, rough, a) * RECIPROCAL_PI;
  float multi  = 0.1159*a;
  return color*(vec3(single) + color*multi);
}

// Normalized Lambert
vec3 DiffuseBRDF(const in IncidentLight directLight, const in GeometricContext geometry, vec3 diffuseColor, float roughnessFactor) {

  vec3 N = geometry.normal;
  vec3 V = geometry.viewDir;
  vec3 L = directLight.direction;
  
  float dotNL = saturate(dot(N,L));
  float dotNV = saturate(dot(N,V));
  vec3 H = normalize(L+V);
  float dotNH = saturate(dot(N,H));
  float dotVH = saturate(dot(V,H));
  float dotLH = saturate(dot(L,H));
  float dotLV = saturate(dot(L,V));
  float a = roughnessFactor * roughnessFactor;
  
  if (diffuseType >= 9) {
    return GGXApproxDiffuse(diffuseColor, dotNL, dotNV, dotLH, dotVH, dotLV, roughnessFactor, a);
  }
  else if (diffuseType >= 8) {
    return GGXOrenNayarDiffuse(diffuseColor, dotNL, dotNV, dotLH, dotVH, dotLV, roughnessFactor, a);
  }
}

vec3 F_Schlick(vec3 specularColor, vec3 H, vec3 V) {
  return (specularColor + (1.0 - specularColor) * pow(1.0 - saturate(dot(V,H)), 5.0));
}

float D_GGX(float a, float dotNH) {
  float a2 = a*a;
  float dotNH2 = dotNH*dotNH;
  float d = dotNH2 * (a2 - 1.0) + 1.0;
  return a2 / (PI * d * d);
}

float G_Smith_Schlick_GGX(float a, float dotNV, float dotNL) {
  float k = a*a*0.5 + EPSILON;
  float gl = dotNL / (dotNL * (1.0 - k) + k);
  float gv = dotNV / (dotNV * (1.0 - k) + k);
  return gl*gv;
}

// Cook-Torrance
vec3 SpecularBRDF(const in IncidentLight directLight, const in GeometricContext geometry, vec3 specularColor, float roughnessFactor) {
  
  vec3 N = geometry.normal;
  vec3 V = geometry.viewDir;
  vec3 L = directLight.direction;
  
  float dotNL = saturate(dot(N,L));
  float dotNV = saturate(dot(N,V));
  vec3 H = normalize(L+V);
  float dotNH = saturate(dot(N,H));
  float dotVH = saturate(dot(V,H));
  float dotLV = saturate(dot(L,V));
  float a = roughnessFactor * roughnessFactor;

  float D = D_GGX(a, dotNH);
  float G = G_Smith_Schlick_GGX(a, dotNV, dotNL);
  vec3 F = F_Schlick(specularColor, V, H);
  return (F*(G*D))/(4.0*dotNL*dotNV+EPSILON);
}

// RenderEquations(RE)
void RE_Direct(const in IncidentLight directLight, const in GeometricContext geometry, const in Material material, inout ReflectedLight reflectedLight) {
  
  float dotNL = saturate(dot(geometry.normal, directLight.direction));
  vec3 irradiance = dotNL * directLight.color;
  
  // punctual light
  irradiance *= PI;
  
  reflectedLight.directDiffuse += irradiance * DiffuseBRDF(directLight, geometry, material.diffuseColor, material.specularRoughness);

  reflectedLight.directSpecular += irradiance * SpecularBRDF(directLight, geometry, material.specularColor, material.specularRoughness);
}

void main() {
  GeometricContext geometry;
  geometry.position = -vViewPosition;
  geometry.normal = normalize(vNormal);
  geometry.viewDir = normalize(vViewPosition);
  
  Material material;
  vec4 temp_dc = texture2D(diffuseMap, vUv);
  vec3 dc = temp_dc.xyz;
  material.diffuseColor = mix(dc, vec3(0.0), metallic);
  material.specularColor = mix(vec3(0.04), dc, metallic);
  vec4 roughnessVector = texture2D(roughnessMap, vUv);
  float roughness = (roughnessVector.x + roughnessVector.y + roughnessVector.z)/3.0;
  material.specularRoughness = roughness/255.0;
  
  // Lighting
  
  ReflectedLight reflectedLight = ReflectedLight(vec3(0.0), vec3(0.0), vec3(0.0), vec3(0.0));
  vec3 emissive = vec3(0.0);
  float opacity = 1.0;
  
  IncidentLight directLight;
  
  // point light
  for (int i=0; i<LIGHT_MAX; ++i) {
    if (i >= numPointLights) break;
    getPointDirectLightIrradiance(pointLights[i], geometry, directLight);
    if (directLight.visible) {
      RE_Direct(directLight, geometry, material, reflectedLight);
    }
  }
  
  // spot light
  for (int i=0; i<LIGHT_MAX; ++i) {
    if (i >= numSpotLights) break;
    getSpotDirectLightIrradiance(spotLights[i], geometry, directLight);
    if (directLight.visible) {
      RE_Direct(directLight, geometry, material, reflectedLight);
    }
  }
  
  // directional light
  for (int i=0; i<LIGHT_MAX; ++i) {
    if (i >= numDirectionalLights) break;
    getDirectionalDirectLightIrradiance(directionalLights[i], geometry, directLight);
    RE_Direct(directLight, geometry, material, reflectedLight);
  }
  
  vec3 outgoingLight = emissive + reflectedLight.directDiffuse + reflectedLight.directSpecular + reflectedLight.indirectDiffuse + reflectedLight.indirectSpecular;
  
  gl_FragColor = vec4(outgoingLight, opacity);
  //gl_FragColor = texture2D(roughnessMap, vUv);

}