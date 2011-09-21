package org.biojava.bio.structure.align.gui;

import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.List;


import javax.swing.Box;
import javax.swing.JButton;

import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JTextField;

import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.ce.CeParameters;
import org.biojava.bio.structure.align.ce.ConfigStrucAligParams;

public class ParameterGUI extends JFrame{

   /**
    * 
    */
   private static final long serialVersionUID = 723386061184110161L;

   ConfigStrucAligParams params ;
   List<Component> textFields;

   
   
   
   @SuppressWarnings("rawtypes")
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
         if ( key.equals(CeParameters.SCORING_STRATEGY) ){
            String[] values = new String[]{"CA only","Sidechain orientation","Angle between sidechains", "CA distance+Angle between sidechains","Sequence Conservation"};
            JComboBox jcbox = new JComboBox(values);
            Integer val = (Integer)value;
            jcbox.setSelectedIndex(val);
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

   @SuppressWarnings({  "rawtypes" })
   protected void setDefault() {
      params.reset();

      List<String> keys  = params.getUserConfigParameters();
      List<Class> types  = params.getUserConfigTypes();
      //List<String> names = params.getUserConfigParameterNames();
      for (int i = 0 ; i < keys.size(); i++){

         Class type = types.get(i);
         Object data = getValue(keys.get(i));
         String name = keys.get(i);
         if ( name.equals(CeParameters.SCORING_STRATEGY)){
            JComboBox field = (JComboBox)  textFields.get(i);
            field.setSelectedIndex((Integer)data);
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
         if ( key.equals(CeParameters.SCORING_STRATEGY)){
            JComboBox field = (JComboBox)  textFields.get(i);
            Integer sel = field.getSelectedIndex();
            value = sel.toString();
         }
         else if ( type == Boolean.class){
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

   @SuppressWarnings({ "unchecked", "rawtypes" })
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
      // first try with get form
      try {
         String methodName = "get" + name;

         @SuppressWarnings("rawtypes")
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
