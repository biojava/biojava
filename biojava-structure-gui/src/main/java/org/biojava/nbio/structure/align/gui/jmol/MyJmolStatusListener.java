/*
 *                    BioJava development code
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
 * Created on Oct 6, 2009
 * Author: Andreas Prlic 
 *
 */

package org.biojava.nbio.structure.align.gui.jmol;

import org.jmol.api.JmolStatusListener;
import org.jmol.constant.EnumCallback;

import javax.swing.*;

import java.util.Hashtable;
import java.util.Map;

public class MyJmolStatusListener implements JmolStatusListener {
	
	private boolean verbose = false;

	JTextField status;
	public MyJmolStatusListener(){
	}
	
	public void setTextField (JTextField statusField) {
		status = statusField;
	}
	
	@Override
	public String createImage(String arg0, String arg1, Object arg2, int arg3) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String eval(String arg0) {
		if(verbose) System.out.println("eval " + arg0);
		return null;
	}

	@Override
	public float[][] functionXY(String arg0, int arg1, int arg2) {
		if(verbose) System.out.println("XY " + arg0 + " " + arg1 + " " + arg2);
		return null;
	}

	@Override
	public float[][][] functionXYZ(String arg0, int arg1, int arg2, int arg3) {
		// TODO Auto-generated method stub
		return null;
	}

	@SuppressWarnings({ "unchecked", "rawtypes" })
	@Override
	public Hashtable getRegistryInfo() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void showUrl(String arg0) {
		status.setText(arg0);

	}

	public void notifyCallback(int arg0, Object[] arg1) {
		status.setText(arg0 +" " + arg1);

	}

	public boolean notifyEnabled(int arg0) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public void setCallbackFunction(String arg0, String arg1) {
		if(verbose) System.out.println("callback:" + arg0 + " " + arg1);
		status.setText(arg0 + " " + arg1);

	}

	public String dialogAsk(String arg0, String arg1) {
		// TODO Auto-generated method stub
		if(verbose) System.out.println("dialogAsk");
		return null;
	}

	public void handlePopupMenu(int arg0, int arg1) {
		// TODO Auto-generated method stub
		if(verbose) System.out.println("handlePopupMenu");
	}

	public void showConsole(boolean arg0) {
		
		// TODO Auto-generated method stub
		if(verbose) System.out.println("showConsole");
		
	}

	
	@Override
	public void notifyCallback(EnumCallback arg0, Object[] arg1) {
		// TODO Auto-generated method stub
		
	}

	
	@Override
	public boolean notifyEnabled(EnumCallback arg0) {
		// TODO Auto-generated method stub
		return false;
	}

	
	@Override
	public Map<String, Object> getProperty(String arg0) {
		// TODO Auto-generated method stub
		return null;
	}

	
		
	public Map<String, Object> getJSpecViewProperty(String arg0) {
		// TODO Auto-generated method stub
		return null;
	}


	@Override
	public void resizeInnerPanel(String data) {
		// TODO Auto-generated method stub
		
	}

	public boolean isVerbose() {
		return verbose;
	}

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

	
	

}
