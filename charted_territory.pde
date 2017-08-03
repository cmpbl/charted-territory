/*

  Author: James Campbell, 2017
  License: MIT
  
  Charted Territory: A simple visualization of geospatial data using Processing.
  
  
  ## KEY VARIABLES
  
  theta -> The positive angle of rotation of the stripes defined by true horizontal (in the easterly direction) and the corresponding rotated line.
           A positive value corresponds to counter-clockwise rotation.
           A zero value would render stripes running in an east-to-west direction; ranks to the north and south of one another.
           
  phi   -> The angle perpendicular to theta.
           phi = theta - PI/2
  
  buffer-> The distance between both stripes and chunk-centers within a stripe.
  
  i     -> The multiplier of buffer which defines the minimum distance to a given stripe (one stripe per value i).
  
  V     -> The view origin -- the point around which the visualization is centered and rotates.       
         
  T     -> The point distance i * buffer from V, along the angle phi.

  C     -> The point at the center of a given stripe-segment or "chunk"
  
  P     -> An evaluated data point.
  
  delta_lat_per_true -> The vertical distance between C and P in the non-rotated coordinate system.
  
  delta_lon_per_true -> The horizontal distance between C and P in the non-rotated coordinate system.
   
  beta  -> The positive angle between the true horizontal (in the easterly direction) and the line CP.
           The final value of beta corresponds to the full positive angle (can be greater than 180).
           Equals atan(delta_lat_per_true / delta_lon_per_true) then modified to equal full positive angle from 0-angle.
  
  alpha -> The positive angle between the rotated horizontal line (in the rotated easterly direction) and the line CP.
           alpha = beta - theta
           
  d     -> The distance between points C and P
           d = sqrt((C_latitude - P_latitude)^2 + (C_longitude - P_longitude)^2)

  delta_lat_per_theta -> The vertical distance between C and P in the coordinate system of the rotated frame of reference.
                         Equals d * sin(alpha)
                         
  delta_lon_per_theta -> The horizontal distance between C and P in the coordinate system of the rotated frame of reference.
                         Equals d * cos(alpha)

*/

// -- CONFIGURATION --

float theta = 0;
int map_offset = 100;
int bg = 252;
int lncolor = 0;

// origin
float origin_lat = 40.718790;
float origin_lon = -73.965345;

// view config
float zoom_factor = 5000;
float vertical_coeff = .2;

float translation_x = 500;
float translation_y = 500;

// stripe parameters
float buffer = .002;
int num_steps = 30;


// -- STRUCTURE DATASET & LOAD --

class Datapoint {
  float lat;
  float lon;
  float value;
  
  Datapoint(float tmpLat, float tmpLon, float tmpValue) {
    lat = tmpLat;
    lon = tmpLon;
    value = tmpValue;
  }
}

ArrayList<Datapoint> loadData() {
  
  ArrayList<Datapoint> data = new ArrayList<Datapoint>();
  
  Table csv = loadTable("dataset.csv", "header");

  for (int i = 0; i < csv.getRowCount(); i++) {
    TableRow row = csv.getRow(i);
    data.add(new Datapoint(row.getFloat("lat"), row.getFloat("lon"), row.getFloat("value")));  
  }
  
  return data;
}

// -- BASEMAP --

void drawMap(){ 
  // Draw the basemap from an appropriately structured CSV form `segment_id|lat|lon`
  
  Table csv = loadTable("nyc_outline.csv", "header");

  stroke(lncolor,100);
  
  int currSeg = -1;
  
  // iterate through coordinate pairs, separating individual segments
  for (int i = 0; i < csv.getRowCount(); i++) {
    TableRow row = csv.getRow(i);
    
    int seg = row.getInt("segment_id");
    float lat = row.getFloat("lat");
    float lon = row.getFloat("lon");
    
    if (currSeg == -1) {
      currSeg = seg; 
      beginShape();
    }
    else if (seg != currSeg) {
      currSeg = seg;
      endShape();
      beginShape();
    }
    
    plotCoord(lat, lon, 0); 
    
  }
  endShape();

}
  
  
// -- TRIG FUNCTIONS --  

float lonToLocalCoord(float lon) {
  // Transform a longitude value to translated and enlarged view.
  
  return (lon-origin_lon)*zoom_factor;
}

float latToLocalCoord(float lat) {
  // Transform a latitude value to translated and enlarged view.
  
  return -(lat-origin_lat)*zoom_factor;
}

void plotCoord(float lat, float lon, float z) {
  // Plot a coordinate and value in translated and enlarged view.
  
  float x = lonToLocalCoord(lon);
  float y = latToLocalCoord(lat);
  vertex(x, y, z * vertical_coeff);
}

void plotLine(float lat1, float lon1, float z1, float lat2, float lon2, float z2) {
  // Plot a line between two coordinate-value tuples in the translated and enlarged view.
  
  float y1 = latToLocalCoord(lat1);
  float x1 = lonToLocalCoord(lon1);
  
  float y2 = latToLocalCoord(lat2);
  float x2 = lonToLocalCoord(lon2);
  
  line (x1, y1, z1, x2, y2, z2);
}

float calculate_lat(float lon, float T_lat, float T_lon, float theta) { 
  // Return the expected latitude value for a given longitude on a specified line of tangency at point T and angle theta.
  
  return tan(theta) * lon + (T_lat - tan(theta) * T_lon);
}


// -- INITIALIZATION -- 
  
void setup()
{
  size(1000, 1000, P3D);
  ortho();
  frameRate(2);
  stroke(lncolor,100);
}

float prior_C_value = 0;
float prior_C_lat = 0;
float prior_C_lon = 0;


// -- MAIN LOOP --

void draw() {   
  // load dataset to chart
  ArrayList<Datapoint> data = loadData();

  // set rotation
  translate(translation_x, translation_y, 0);
  //rotateY(radians(PI/10));
  rotateX(radians(70));
  rotateZ(theta-PI/10);
  theta+=PI*.005;
  if (theta > 2*PI) {
    theta = theta-2*PI;
  }
  
  // clear between frames and redraw base map
  background(bg);
  noFill();
  drawMap();
 
  // initialize variables
  float lat;
  float lon;
  float value;
  float T_lat;
  float T_lon;
  float phi;
  
  // set phi perpendicular to theta
  phi = theta + PI/2;
  
  // iterate through step-count below and above origin to determine stripe membership
  for (int i = -num_steps; i < num_steps; i++) {
  
    // calculate T point of tangency
    T_lat = origin_lat + sin(phi) * i * buffer;
    T_lon = origin_lon + cos(phi) * i * buffer;
    
    // calculate incremental offsets for latitude and longitude
    float shift_step = buffer * abs(sin(phi));
    float apparent_depth = buffer * abs(cos(phi));

    // initialize stripe object to store points along current step's tangent line
    ArrayList<Datapoint> stripe = new ArrayList<Datapoint>();
    
    // evaluate each point in dataset for membership in stripe
    for (int m = 0; m < data.size(); m++) {
      lat = data.get(m).lat;
      lon = data.get(m).lon;
      value = data.get(m).value;
      
      float predicted_lat = calculate_lat(lon, T_lat, T_lon, theta);
      float diff = abs(lat - predicted_lat);
      
      // add qualifying datapoints to stripe
      if (diff < buffer) { 
        stripe.add(new Datapoint(lat, lon, value));
      }
    }
    
    // iterate through step-count left and right of origin to determine chunk membership
    for (int k = -num_steps; k < num_steps; k++) {
  
      // initialize chunk object
      ArrayList<Datapoint> chunk = new ArrayList<Datapoint>();
      
      // calculate center point of the current chunk
      float C_lon = T_lon + (k * shift_step);
      float C_lat = calculate_lat(C_lon, T_lat, T_lon, theta);  
  
      // reset ridge line for new chunk or discontinuity
      if(prior_C_lat == 0) {
        prior_C_lat = C_lat;
        prior_C_lon = C_lon;
        prior_C_value = 0;  
      }
      
      // evaluate each stripe-point for membership in chunk
      for (int j = 0; j < stripe.size(); j++){
        // compare longitude to determine if in chunk -- i already know it's in the stripe
        lat = stripe.get(j).lat;
        lon = stripe.get(j).lon;
        value = stripe.get(j).value;
        
        float delta_lat_per_true = C_lat - lat;
        float delta_lon_per_true = C_lon - lon;
        
        float d = sqrt(pow(delta_lat_per_true,2) + pow(delta_lon_per_true,2));
        
        float small_beta = atan(delta_lat_per_true / delta_lon_per_true);
        float beta = 0;
        
        // need to catch ='s
        if (delta_lat_per_true > 0 && delta_lon_per_true > 0) {  // quadrant 1 (unrotated)
          beta = small_beta;
        }
         else if (delta_lat_per_true > 0 && delta_lon_per_true < 0) { // quadrant 2
          beta = PI - small_beta;
        }
         else if (delta_lat_per_true < 0 && delta_lon_per_true < 0) { // quadrant 3
          beta = PI + small_beta;
        }
        else { // quadrant 4
          beta = 2*PI - small_beta;
        }
        
        float alpha = beta - theta;
        float delta_lat_per_theta = d * sin(alpha);
        float delta_lon_per_theta = d * cos(alpha);
        
        // add qualifying datapoints to chunk
        if (abs(delta_lat_per_theta) < buffer/2 && abs(delta_lon_per_theta) < buffer/2) {
          chunk.add(new Datapoint(lat, lon, value));
        }
      }
      
      // aggregate chunk values
      value = 0;
      for (int j = 0; j < chunk.size(); j++) {
         value += chunk.get(j).value;
      }
      
      // plot resulting coordinate-value for chunk
      if (prior_C_value > 0 || value > 0) {
        fill(252);
        beginShape();
        plotCoord(prior_C_lat, prior_C_lon, 0); 
        plotCoord(prior_C_lat, prior_C_lon, prior_C_value); 
        plotCoord(C_lat, C_lon, value);
        plotCoord(C_lat, C_lon, 0);
        endShape();
      }
      
      // set base for next ridge iteration
      prior_C_lat = C_lat;
      prior_C_lon = C_lon;
      prior_C_value = value;
    } 
    
    // reset for a new ridge
    prior_C_lat = 0;
  }

}