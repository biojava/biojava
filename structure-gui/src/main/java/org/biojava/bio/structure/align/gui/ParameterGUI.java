package org.biojava.bio.structure.align.gui;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.List;


import javax.swing.Box;
import javax.swing.JButton;

import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JTextField;

import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.ce.ConfigStrucAligParams;

public class ParameterGUI extends JFrame{

	/**
	 * 
	 */
	private static final long serialVersionUID = 723386061184110161L;

	ConfigStrucAligParams params ;
	List<JTextField> textFields;
	
	public ParameterGUI(StructureAlignment alignment){

		ConfigStrucAligParams params = alignment.getParameters();
		
		if ( params == null)
			return;
		this.params = params;
		
		String method = alignment.getAlgorithmName();
		this.setTitle("Parameters for " + method);
		
		
		List<String> names = params.getUserConfigParameterNames();
		List<String> keys  = params.getUserConfigParameters();
		List<Class> types  = params.getUserConfigTypes();
		
		List<String> helps = params.getUserConfigHelp();
		textFields = new ArrayList<JTextField>();
		Box vBox = Box.createVerticalBox();

		for (int i = 0 ; i < names.size(); i++){
		    Class type = types.get(i);
		    
			Box hBox = Box.createHorizontalBox();

			JLabel label = new JLabel(names.get(i));
			String help = helps.get(i);
			label.setToolTipText(help);

			Object value = getValue(keys.get(i));

			String data = value.toString();
			JTextField field = new JTextField(10);
			if ( type == String[].class) {
			   String stuff = "";
               for ( String da : (String[]) value){
                  stuff += da + " ";
               }
               data = stuff;
              
               
			}
			field.setText(data);
			field.setToolTipText(help);

			hBox.add(label);
			hBox.add(Box.createGlue());
			hBox.add(field);

			vBox.add(hBox);
			
			textFields.add(field);

		}


		JButton abort = new JButton("Cancel");
		abort.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent event) {
				destroy();
				dispose();	         }
		});
		
		JButton defaultB = new JButton("Default");
		defaultB.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent event) {
				setDefault();
			}
		});

		JButton close = new JButton("Apply");

		close.addActionListener(new ActionListener(){
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
	
	protected void setDefault() {
		params.reset();
		
		List<String> keys = params.getUserConfigParameters();
		List<Class> types = params.getUserConfigTypes();
		for (int i = 0 ; i < keys.size(); i++){
			JTextField field = textFields.get(i);
			Class type = types.get(i);
			Object data = getValue(keys.get(i));
			if ( type.isArray()){
			   String stuff = "";
			   for ( String da : (String[]) data){
			      stuff += da + " ";
			   }
			   System.out.println(type + "setting string array:" + stuff);
			   field.setText(stuff);
			} else {
			   System.out.println(type + "setting string array:" + data.toString());
			   field.setText(data.toString()); 
			}
			
			field.updateUI();
		}
		this.repaint();
		
	}

	private void destroy(){
		//avoid memory leaks...
		textFields = null;
		params = null;
	}

	@SuppressWarnings("unchecked")
	protected void storeParameters() {
		List<String> names = params.getUserConfigParameterNames();
		List<String> keys = params.getUserConfigParameters();
		List<Class> types = params.getUserConfigTypes();
		
		for (int i = 0 ; i < names.size(); i++){
			JTextField field = textFields.get(i);
			String value = field.getText();
			
			Class type = types.get(i);
			String key  = keys.get(i);
			setValue(key, type, value);
		}
		
		System.out.println("new parameters: " + params.toString());

	}

	@SuppressWarnings("unchecked")
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
			} else if ( type == String[].class) {
			   data = value.split(" ");
			   
			}
			
			if (data == null){
				System.err.println("Could not set value " + value + " for field " + name);
				return;
			}
			 m.invoke(params, data);

			
		} catch (Exception e){
			e.printStackTrace();
		
		}
		
	}

	@SuppressWarnings("unchecked")
	private Object  getValue(String name){

		try {
			String methodName = "get" + name;

			Class paramC = params.getClass();

			Method m =paramC.getMethod(methodName,null);

			Object value = m.invoke(params);

			return value;
		} catch (Exception e){
			e.printStackTrace();
			return null;
		}


	}
}
