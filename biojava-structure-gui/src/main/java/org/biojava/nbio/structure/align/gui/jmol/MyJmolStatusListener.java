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

import java.util.Map;

import javax.swing.JTextField;

import org.jmol.api.JmolStatusListener;
import org.jmol.c.CBK;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class MyJmolStatusListener implements JmolStatusListener {
	private static final Logger logger = LoggerFactory.getLogger(MyJmolStatusListener.class);

	JTextField status;
	public MyJmolStatusListener(){
	}

	public void setTextField (JTextField statusField) {
		status = statusField;
	}

	@Override
	public String createImage(String arg0, String arg1, Object arg2, int arg3) {
		return null; //Cancelled
	}

	@Override
	public String eval(String arg0) {
		logger.debug("eval {}",arg0);
		return null;
	}

	@Override
	public float[][] functionXY(String arg0, int arg1, int arg2) {
		logger.debug("XY {} {} {}",arg0,arg1, arg2);
		return null; //Ignore isosurface commands
	}

	@Override
	public float[][][] functionXYZ(String arg0, int arg1, int arg2, int arg3) {
		logger.debug("XYZ {} {} {} {}",arg0,arg1, arg2, arg3);
		return null; //Ignore isosurface commands
	}

	@Override
	public Map<String, Object> getRegistryInfo() {
		return null; //Ignore
	}

	@Override
	public void showUrl(String arg0) {
		status.setText(arg0);

	}

	public void notifyCallback(int arg0, Object[] arg1) {
		status.setText(arg0 +" " + arg1);

	}

	public boolean notifyEnabled(int arg0) {
		return false;
	}

	@Override
	public void setCallbackFunction(String arg0, String arg1) {
		logger.debug("callback: {} {}", arg0, arg1);
		status.setText(arg0 + " " + arg1);
	}

	public String dialogAsk(String arg0, String arg1) {
		logger.debug("dialogAsk {} {}",arg0,arg1);
		return null; //Ignore
	}

	public void handlePopupMenu(int arg0, int arg1) {
		logger.debug("handlePopupMenu {} {}",arg0,arg1);
	}

	public void showConsole(boolean arg0) {
		logger.debug("showConsole {}",arg0);
	}



	@Override
	public void notifyCallback(CBK message, Object[] data) {
	}

	@Override
	public boolean notifyEnabled(CBK type) {
		return false;
	}

	@Override
	public Map<String, Object> getJSpecViewProperty(String arg0) {
		return null;
	}


	@Override
	public int[] resizeInnerPanel(String data) {
		return null;
	}

}
