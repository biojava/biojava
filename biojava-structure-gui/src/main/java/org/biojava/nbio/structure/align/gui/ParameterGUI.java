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
 */
package org.biojava.nbio.structure.align.gui;

import org.biojava.nbio.structure.align.ce.ConfigStrucAligParams;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.List;

/**
 * UI for {@link ConfigStrucAligParams}, for the AlignmentGUI.
 *
 * Visually parameters are displayed as a label and an input field such as a text box.
 * The following methods are used to determine which properties are accessible:
 * <ul>
 * <li> {@link ConfigStrucAligParams#getUserConfigParameterNames() getUserConfigParameterNames()}
 *      Parameter labels
 * <li> {@link ConfigStrucAligParams#getUserConfigParameters() getUserConfigParameters()}
 *      Parameter method names, to be prepended with 'set' or 'get'
 * <li> {@link ConfigStrucAligParams#getUserConfigTypes() getUserConfigTypes()}
 *      Types for each parameter (Integer, Float, Boolean, Enum, String[])
 * <li> {@link ConfigStrucAligParams#getUserConfigHelp() getUserConfigHelp()}
 *      Alt text for each parameter
 * </ul>
 *
 * @author Andreas Prlic
 *
 */
public class ParameterGUI extends JFrame{

	private static final long serialVersionUID = 723386061184110161L;

	private ConfigStrucAligParams params ;
	private List<Component> textFields;

	/**
	 * Constructor for a ParameterGUI. Generalization for any type of
	 * Structural Alignment algorithm that implements the parameter interface.
	 *
	 * @param params parameter bean
	 * @param algorithm name of the algorithm
	 */
	@SuppressWarnings("rawtypes")
	public ParameterGUI(ConfigStrucAligParams params, String algorithm) {

		if (params == null) return;
		this.params = params;

		this.setTitle("Parameters for " + algorithm);

		List<String> names = params.getUserConfigParameterNames();
		List<String> keys  = params.getUserConfigParameters();
		List<Class> types  = params.getUserConfigTypes();

		List<String> helps = params.getUserConfigHelp();

		// quick check for bugs in params
		assert(names.size() == keys.size());
		assert(names.size() == types.size());
		assert(names.size() == helps.size());

		textFields = new ArrayList<Component>();
		Box vBox = Box.createVerticalBox();

		for (int i = 0 ; i < keys.size(); i++){
			Class type = types.get(i);

			Box hBox = Box.createHorizontalBox();
			String name = names.get(i);
			JLabel label = new JLabel(name);
			String help = helps.get(i);
			label.setToolTipText(help);
			String key = keys.get(i);
			Object value = getValue(key);

			String data = value.toString();
			Component field;
			if ( type.isEnum() ) {
				Object[] values = type.getEnumConstants();
				JComboBox jcbox = new JComboBox(values);
				jcbox.setSelectedItem(value);
				field = jcbox;

			} else if ( type == Boolean.class){

				String[] values = new String[]{"true","false"};
				JComboBox jcbox = new JComboBox(values);
				if ( data.equalsIgnoreCase("false"))
					jcbox.setSelectedIndex(1);
				else
					jcbox.setSelectedIndex(0);

				field = jcbox;

				//field.setToolTipText(help);

			} else {
				JTextField tfield = new JTextField(10);

				if ( type == String[].class) {
					String stuff = "";
					for ( String da : (String[]) value){
						stuff += da + " ";
					}
					data = stuff;


				}
				tfield.setText(data);
				tfield.setToolTipText(help);
				field = tfield;
			}

			hBox.add(label);
			hBox.add(Box.createGlue());
			hBox.add(field);

			vBox.add(hBox);

			textFields.add(field);

		}


		JButton abort = new JButton("Cancel");
		abort.addActionListener(new ActionListener(){
			@Override
			public void actionPerformed(ActionEvent event) {
				destroy();
				dispose();	         }
		});

		JButton defaultB = new JButton("Default");
		defaultB.addActionListener(new ActionListener(){
			@Override
			public void actionPerformed(ActionEvent event) {
				setDefault();
			}
		});

		JButton close = new JButton("Apply");

		close.addActionListener(new ActionListener(){
			@Override
			public void actionPerformed(ActionEvent event) {

				storeParameters();

				destroy();
				dispose();	         }
		});

		Box hBox = Box.createHorizontalBox();
		hBox.add(abort);
		hBox.add(Box.createGlue());
		hBox.add(defaultB);
		hBox.add(Box.createGlue());
		hBox.add(close);

		vBox.add(hBox);
		this.getContentPane().add(vBox);
		this.pack();
		this.setVisible(true);
	}

	@SuppressWarnings({  "rawtypes" })
	protected void setDefault() {
		params.reset();

		List<String> keys  = params.getUserConfigParameters();
		List<Class> types  = params.getUserConfigTypes();
		//List<String> names = params.getUserConfigParameterNames();
		for (int i = 0 ; i < keys.size(); i++){

			Class type = types.get(i);
			Object data = getValue(keys.get(i));
			if( type.isEnum()) {
				JComboBox field = (JComboBox)  textFields.get(i);
				field.setSelectedItem(data);
				field.updateUI();
			}  else if ( type == Boolean.class){
				JComboBox field = (JComboBox)  textFields.get(i);
				if ( data.toString().equalsIgnoreCase("false"))
					field.setSelectedIndex(1);
				else
					field.setSelectedIndex(0);
				field.updateUI();

			} else {
				JTextField field = (JTextField)textFields.get(i);
				if ( type.isArray()){
					String stuff = "";
					for ( String da : (String[]) data){
						stuff += da + " ";
					}

					field.setText(stuff);
				} else {

					field.setText(data.toString());
				}
				field.updateUI();
			}


		}
		this.repaint();

	}

	private void destroy(){
		//avoid memory leaks...
		textFields = null;
		params = null;
	}

	@SuppressWarnings("rawtypes")
	protected void storeParameters() {
		//List<String> names = params.getUserConfigParameterNames();
		List<String> keys = params.getUserConfigParameters();
		List<Class> types = params.getUserConfigTypes();

		for (int i = 0 ; i < keys.size(); i++){
			Class type = types.get(i);
			String key  = keys.get(i);
			// String name = keys.get(i);
			String value = null;
			System.out.println(key);
			if( type.isEnum() ) {
				JComboBox field = (JComboBox)  textFields.get(i);
				Enum sel = (Enum)field.getSelectedItem();
				value = sel.name();
			} else if ( type == Boolean.class){
				JComboBox field = (JComboBox)  textFields.get(i);
				int sel = field.getSelectedIndex();
				Boolean flag = true;
				if ( sel == 1 )
					flag = false;
				value = flag.toString();
			} else {
				JTextField field = (JTextField)textFields.get(i);
				value = field.getText();
			}

			setValue(key, type, value);
		}

		System.out.println("new parameters: " + params.toString());

	}

	@SuppressWarnings({ "unchecked" })
	private void setValue(String name, Class type, String value) {
		try {
			String methodName = "set" + name;

			Class paramC = params.getClass();

			Method m =paramC.getMethod(methodName,type);


			Object data = null;

			if ( type == Integer.class){
				data = Integer.parseInt(value);
			} else if ( type == Double.class){
				data = Double.parseDouble(value);
			} else if ( type == Float.class) {
				data = Float.parseFloat(value);
			} else if ( type == Boolean.class) {
				data = Boolean.parseBoolean(value);
			} else if ( type == Short.class) {
				data = Short.parseShort(value);
			} else if ( type == String[].class) {
				data = value.split(" ");
			} else if ( type.isEnum() ) {
				data = Enum.valueOf(type, value);
			}

			if (data == null){
				System.err.println("Could not set value " + value +
						" for field " + name);
				return;
			}
			m.invoke(params, data);


		} catch (Exception e){
			e.printStackTrace();

		}

	}

	@SuppressWarnings("unchecked")
	private Object  getValue(String name){
		// first try with get form
		try {
			String methodName = "get" + name;

			Class paramC = params.getClass();

			Method m;
			try {
				//try boolean getter
				m = paramC.getMethod(methodName,(Class[])null);
			} catch(NoSuchMethodException e) {
				//try boolean getter
				methodName = "is" + name;
				m = paramC.getMethod(methodName,(Class[])null);
			}

			Object value = m.invoke(params);

			return value;
		} catch (Exception e){
			e.printStackTrace();
			return null;
		}


	}
}
