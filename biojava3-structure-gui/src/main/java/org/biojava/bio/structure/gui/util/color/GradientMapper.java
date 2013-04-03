/*
 *                  BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 * 
 * Created on Aug 3, 2007
 */
package org.biojava.bio.structure.gui.util.color;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.color.ColorSpace;
import java.util.Collection;
import java.util.Map;
import java.util.NavigableMap;
import java.util.Set;
import java.util.TreeMap;

import javax.swing.JFrame;
import javax.swing.JPanel;

import org.biojava.bio.structure.gui.util.color.LinearColorInterpolator.InterpolationDirection;

/**
 * Maps a set of real values onto a gradient.
 * 
 * The real line is partitioned into segments [a,b). The endpoint of each segment is labeled with a color.
 * Colors are linearly interpolated between finite endpoints. Endpoints implicitly exist for 
 * Double.NEGATIVE_INFINITY and Double.POSITIVE_INFINITY, representing default colors. Thus any point
 * in the segment [-Inf,a) is labeled with the negInf color, and any point in [b,Inf] is labeled with the posInf
 * color. If no endpoints are present, the posInf color is used as default.
 * 
 * Common gradients are predefined an may be instantiated through 
 * GradientMapper.getGradientMapper().
 * 
 * @author Spencer Bliven
 *
 */
public class GradientMapper implements ContinuousColorMapper, Map<Double, Color> {

	public static final int BLACK_WHITE_GRADIENT = 1;
	public static final int WHITE_BLACK_GRADIENT = 2;
	public static final int RED_BLUE_GRADIENT = 3;
	public static final int RAINBOW_GRADIENT = 4;
	public static final int RAINBOW_INTENSITY_GRADIENT = 5;
	
	private NavigableMap<Double,Color> mapping;
	private ColorInterpolator interpolator;
	
	public GradientMapper() {
		this(Color.black,Color.white);
	}
	public GradientMapper(Color negInf, Color posInf) {
		this(negInf,posInf,ColorSpace.getInstance(ColorSpace.CS_sRGB));
	}
	public GradientMapper(Color negInf, Color posInf, ColorSpace cspace) {
		mapping = new TreeMap<Double,Color>();
		mapping.put(Double.NEGATIVE_INFINITY, negInf);
		mapping.put(Double.POSITIVE_INFINITY, posInf);
		interpolator = new LinearColorInterpolator(cspace);
	}

	/**
	 * Constructs a gradientMapper to draw one of the pre-defined gradients
	 * 
	 * For example,
	 * GradientMapper.getGradientMapper(GradientMapper.RAINBOW_GRADIENT, 0, 10)
	 * 
	 * @param gradientType One of the gradient types, eg GradientMapper.BLACK_WHITE_GRADIENT
	 * @param min Start of the gradient
	 * @param max End of the gradient
	 * @return
	 */
	public static GradientMapper getGradientMapper(int gradientType, double min, double max) {
		GradientMapper gm;
		switch( gradientType ) {
		case BLACK_WHITE_GRADIENT:
			gm = new GradientMapper(Color.BLACK, Color.WHITE);
			gm.put(min, Color.BLACK);
			gm.put(max, Color.WHITE);
			return gm;
		case WHITE_BLACK_GRADIENT:
			gm = new GradientMapper(Color.WHITE, Color.BLACK);
			gm.put(min, Color.WHITE);
			gm.put(max, Color.BLACK);
			return gm;
		case RED_BLUE_GRADIENT:
			gm = new GradientMapper(Color.RED, Color.BLUE);
			gm.put(min, Color.RED);
			gm.put(max, Color.BLUE);
			return gm;
		case RAINBOW_GRADIENT: {
			//Set up interpolation in HSV colorspace
			ColorSpace hsv = HSVColorSpace.getHSVColorSpace();
			LinearColorInterpolator interp = new LinearColorInterpolator(hsv);
			interp.setInterpolationDirection(0, InterpolationDirection.UPPER);
			
			Color hsvLow = new Color(hsv,new float[] {0f, 1f, 1f},1f);
			Color hsvHigh = new Color(hsv,new float[] {1f, 1f, 1f},1f);
			
			gm = new GradientMapper(hsvLow, hsvHigh, hsv);
			gm.put(min, hsvLow);
			gm.put(max, hsvHigh);
			gm.setInterpolator(interp);
			return gm;
		}
		case RAINBOW_INTENSITY_GRADIENT: {
			//Set up interpolation in HSV colorspace
			ColorSpace hsv = HSVColorSpace.getHSVColorSpace();
			LinearColorInterpolator interp = new LinearColorInterpolator(hsv);
			interp.setInterpolationDirection(0, InterpolationDirection.LOWER);
			
			Color hsvLow = new Color(hsv,new float[] {1f, 1f, 1f},1f);
			Color hsvHigh = new Color(hsv,new float[] {0f, 1f, 0f},1f);
			
			gm = new GradientMapper(hsvLow, hsvHigh, hsv);
			gm.put(min, hsvLow);
			gm.put(max, hsvHigh);
			gm.setInterpolator(interp);
			return gm;
		}
		default:
			throw new IllegalArgumentException("Unsupported gradient "+gradientType);
		}
	}
	/**
	 * @param value
	 * @return
	 * @see org.biojava.bio.structure.gui.util.color.ContinuousColorMapper#getColor(double)
	 */
	public Color getColor(double value) {
		Double left = mapping.floorKey(value);
		Double right = mapping.higherKey(value);
		
		//don't interpolate to infinity
		if(right == null || right.isInfinite()) {
			return mapping.get(Double.POSITIVE_INFINITY);
		}
		if(left == null || left.isInfinite()) {
			return mapping.get(Double.NEGATIVE_INFINITY);
		}
		
		// fraction of left color to use
		float alpha = (float) ((right-value)/(right-left));
		return interpolator.interpolate(mapping.get(left),mapping.get(right),alpha);
	}

	
	
	/*-*************************
	 *  Map methods
	 ***************************/


	/**
	 * Clears all finite endpoints
	 * 
	 * @see java.util.Map#clear()
	 */
	public void clear() {
		Color neg = mapping.get(Double.NEGATIVE_INFINITY);
		Color pos = mapping.get(Double.POSITIVE_INFINITY);
		mapping.clear();		
		mapping.put(Double.NEGATIVE_INFINITY, neg);
		mapping.put(Double.POSITIVE_INFINITY, pos);
	}
	/**
	 * @param position
	 * @return
	 * @see java.util.Map#containsKey(java.lang.Object)
	 */
	public boolean containsKey(Object position) {
		return mapping.containsKey(position);
	}
	/**
	 * @param color
	 * @return
	 * @see java.util.Map#containsValue(java.lang.Object)
	 */
	public boolean containsValue(Object color) {
		return mapping.containsValue(color);
	}
	/**
	 * @return
	 * @see java.util.Map#entrySet()
	 */
	public Set<java.util.Map.Entry<Double, Color>> entrySet() {
		return mapping.entrySet();
	}
	/**
	 * @param position
	 * @return The color of the endpoint at position, or null if no endpoint exists there
	 * @see java.util.Map#get(java.lang.Object)
	 */
	public Color get(Object position) {
		return mapping.get(position);
	}
	/**
	 * @return true if this gradient does not contain finite endpoints
	 * @see java.util.Map#isEmpty()
	 */
	public boolean isEmpty() {
		return mapping.size() <= 2;
	}
	/**
	 * @return
	 * @see java.util.Map#keySet()
	 */
	public Set<Double> keySet() {
		return mapping.keySet();
	}
	/**
	 * Adds a gradient endpoint at the specified position.
	 * @param position The endpoint position. May be Double.POSITIVE_INFINITY or Double.NEGATIVE_INFINITY for endpoints.
	 * @param color
	 * @return
	 * @see java.util.Map#put(java.lang.Object, java.lang.Object)
	 */
	public Color put(Double position, Color color) {
		if( position == null ) {
			throw new NullPointerException("Null endpoint position");
		}
		if( color == null ){
			throw new NullPointerException("Null colors are not allowed.");
		}
		return mapping.put(position, color);
	}
	/**
	 * @param m
	 * @see java.util.Map#putAll(java.util.Map)
	 */
	public void putAll(Map<? extends Double, ? extends Color> m) {
		mapping.putAll(m);		
	}
	/**
	 * @param position
	 * @return
	 * @see java.util.Map#remove(java.lang.Object)
	 */
	public Color remove(Object position) {
		if( ((Double)position).isInfinite() ) {
			throw new UnsupportedOperationException("Cannot remove infinite endpoints");
		}
		return mapping.remove(position);
	}
	/**
	 * @return Number of finite endpoints
	 * @see java.util.Map#size()
	 */
	public int size() {
		return mapping.size()-2;
	}
	/**
	 * @return
	 * @see java.util.Map#values()
	 */
	public Collection<Color> values() {
		return mapping.values();
	}
	
	
	/**
	 * @return the interpolator
	 */
	public ColorInterpolator getInterpolator() {
		return interpolator;
	}
	/**
	 * @param interpolator the interpolator to set
	 */
	public void setInterpolator(ColorInterpolator interpolator) {
		this.interpolator = interpolator;
	}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		GradientMapper[] mappers = new GradientMapper[20];
		int i = 0;
		ColorSpace hsv = HSVColorSpace.getHSVColorSpace();
		LinearColorInterpolator interp;
		
		
		// RGB colorspace
		mappers[i] = new GradientMapper(Color.black, Color.white);
		mappers[i].put(-5., Color.red);
		mappers[i].put(5., Color.blue);
		i++;
		
		// Premade
		mappers[i] = GradientMapper.getGradientMapper(BLACK_WHITE_GRADIENT,-5,5);
		i++;
		mappers[i] = GradientMapper.getGradientMapper(RAINBOW_INTENSITY_GRADIENT,-5,5);
		//mappers[i].put(Double.NEGATIVE_INFINITY, mappers[i].get(Double.NEGATIVE_INFINITY).brighter());
		//mappers[i].put(Double.POSITIVE_INFINITY, mappers[i].get(Double.POSITIVE_INFINITY).darker());
		i++;
		
		// Rainbow
		mappers[i] = new GradientMapper(Color.black, Color.white, hsv);
		mappers[i].put(-5., new Color(hsv,new float[] {0f, 1f, 1f},1f));
		mappers[i].put( 5., new Color(hsv,new float[] {1f, 1f, 1f},1f));
		i++;

		// HSV INNER
		mappers[i] = new GradientMapper(Color.black, Color.white, hsv);
		mappers[i].put( 5., Color.red);
		mappers[i].put(-5., Color.blue);
		i++;
		
		// HSV OUTER
		interp = new LinearColorInterpolator(hsv);
		interp.setInterpolationDirection(0, InterpolationDirection.OUTER);
		mappers[i] = new GradientMapper(Color.black, Color.white, hsv);
		mappers[i].put( 5., Color.red);
		mappers[i].put(-5., Color.blue);
		mappers[i].setInterpolator(interp);
		i++;
		
		// HSV UPPER
		interp = new LinearColorInterpolator(hsv);
		interp.setInterpolationDirection(0, InterpolationDirection.UPPER);
		mappers[i] = new GradientMapper(Color.black, Color.white, hsv);
		mappers[i].put( 5., Color.red);
		mappers[i].put(-5., Color.blue);
		mappers[i].setInterpolator(interp);
		i++;
		
		// HSV LOWER
		interp = new LinearColorInterpolator(hsv);
		interp.setInterpolationDirection(0, InterpolationDirection.LOWER);
		mappers[i] = new GradientMapper(Color.black, Color.white, hsv);
		mappers[i].put( 5., Color.red);
		mappers[i].put(-5., Color.blue);
		mappers[i].setInterpolator(interp);
		i++;

		// Mimic DefaultMapper
		interp = new LinearColorInterpolator(hsv);
		interp.setInterpolationDirection(0, InterpolationDirection.INNER);
		mappers[i] = new GradientMapper(Color.green, Color.black, hsv);
		mappers[i].put( 0., new Color(hsv,new float[] {1f, .9f, 1f},1f));
		mappers[i].put(10., new Color(hsv,new float[] {0f, .9f, 0f},1f));
		mappers[i].setInterpolator(interp);
		i++;
		
		// Better DefaultGradient
		interp = new LinearColorInterpolator(hsv);
		interp.setInterpolationDirection(0, InterpolationDirection.INNER);
		mappers[i] = new GradientMapper(Color.green, Color.black, hsv);
		mappers[i].put( 0., new Color(hsv,new float[] {1f, .9f, 1f},1f));
		mappers[i].put( 1., new Color(hsv,new float[] {0f, .9f, 1f},1f));
		mappers[i].put( 1+1e-6, Color.white);
		mappers[i].put(10., Color.black);
		mappers[i].setInterpolator(interp);
		i++;
		// Better DefaultGradient
		interp = new LinearColorInterpolator(hsv);
		interp.setInterpolationDirection(0, InterpolationDirection.INNER);
		mappers[i] = new GradientMapper(Color.green, Color.black, hsv);
		mappers[i].put( 0., new Color(hsv,new float[] {1f, .9f, 1f},1f));
		mappers[i].put( 1., new Color(hsv,new float[] {.2f, .9f, 1f},1f));
		mappers[i].put( 1+1e-6, Color.white);
		mappers[i].put(10., Color.black);
		mappers[i].setInterpolator(interp);
		i++;
		
		
		DefaultMatrixMapper defaultMapper = new DefaultMatrixMapper(10f,.9f);
		
		
		
		JFrame frame = new JFrame("GradientMapper");
		JPanel main = new JPanel();
		main.setPreferredSize(new Dimension(300,500));

		for(int j=0;j<i;j++) {
			GradientPanel grad1 = new GradientPanel(mappers[j],-10,10);
			//grad1.setPreferredSize(new Dimension(500,50));
			main.add(grad1);
		}
		GradientPanel grad2 = new GradientPanel(defaultMapper,-10,10);
		//grad2.setPreferredSize(new Dimension(500,50));
		main.add(grad2);
		//main.add(new GradientPanel(defaultMapper,-10,10));


		frame.getContentPane().add(main);
		frame.pack();
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.setVisible(true);
		
	}


}
