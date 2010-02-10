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

package org.biojava.bio;

import java.util.ArrayList;
import java.util.List;

import junit.framework.TestCase;

import org.biojava.bio.symbol.Location;

/**
 * Tests for AnnotationType.
 *
 * @author Matthew Pocock
 * @author Thomas Down
 * @since 1.3
 */
public class AnnotationTypeTest
extends TestCase {
  protected Location cc_zero;
  protected Location cc_zero_or_one;
  protected Location cc_one;
  protected Location cc_one_or_more;
  protected Location cc_any;

  protected PropertyConstraint pc_any;
  protected PropertyConstraint pc_class_string;
  protected PropertyConstraint pc_class_list;
  protected PropertyConstraint pc_class_arraylist;
  
  protected PropertyConstraint pc_exact_fish;
  protected PropertyConstraint pc_exact_dog;
  
  protected PropertyConstraint pc_enum_animals;
  protected PropertyConstraint pc_enum_birds;
  protected PropertyConstraint pc_enum_primary_colors;
  protected PropertyConstraint pc_enum_warm_colors;
  
  protected AnnotationType at_any;
  protected AnnotationType at_color;
  protected AnnotationType at_colors;
  protected AnnotationType at_colorsish;
  protected AnnotationType at_color_name;
  protected AnnotationType at_color_dog;
  
  public AnnotationTypeTest(String name) {
    super(name);
  }
  
  protected void setUp() throws Exception {
    cc_zero = CardinalityConstraint.ZERO;
    cc_zero_or_one = CardinalityConstraint.ZERO_OR_ONE;
    cc_one = CardinalityConstraint.ONE;
    cc_one_or_more = CardinalityConstraint.ONE_OR_MORE;
    cc_any = CardinalityConstraint.ANY;
    
    pc_any = PropertyConstraint.ANY;
    pc_class_string = new PropertyConstraint.ByClass(String.class);
    pc_class_list = new PropertyConstraint.ByClass(List.class);
    pc_class_arraylist = new PropertyConstraint.ByClass(ArrayList.class);
    
    pc_exact_fish = new PropertyConstraint.ExactValue("fish");
    pc_exact_dog = new PropertyConstraint.ExactValue("dog");
    
    pc_enum_animals = new PropertyConstraint.Enumeration(new Object[] { "dog", "cat", "antelope", "crow", "hawk", "tit" });
    pc_enum_birds = new PropertyConstraint.Enumeration(new Object[] { "crow", "hawk", "tit" });
    pc_enum_primary_colors = new PropertyConstraint.Enumeration(new Object[] { "red", "green", "blue" });
    pc_enum_warm_colors = new PropertyConstraint.Enumeration(new Object[] {"red", "orange", "russet"});
    
    at_any = AnnotationType.ANY;
    
    at_color = new AnnotationType.Impl();
    at_color.setConstraints("color", pc_enum_primary_colors, cc_one);

    at_colors = new AnnotationType.Impl();
    at_colors.setConstraints("color", pc_enum_primary_colors, cc_one_or_more);
    
    at_colorsish = new AnnotationType.Impl();
    at_colorsish.setConstraints("color", pc_enum_primary_colors, cc_any);
    
    at_color_name = new AnnotationType.Impl();
    at_color_name.setConstraints("color", pc_enum_primary_colors, cc_one);
    at_color_name.setConstraints("name", pc_class_string, cc_one);
    
    at_color_dog = new AnnotationType.Impl();
    at_color_dog.setConstraints("color", pc_enum_primary_colors, cc_one);
    at_color_dog.setConstraints("name", pc_exact_dog, cc_one);
  }
  
  public void testClassProperties() {
    // any
    assertTrue("are subconstraint: " + pc_any + ", " + pc_any, pc_any.subConstraintOf(pc_any));
    assertTrue("are subconstraint: " + pc_any + ", " + pc_class_string, pc_any.subConstraintOf(pc_class_string));
    assertTrue("are subconstraint: " + pc_any + ", " + pc_class_list, pc_any.subConstraintOf(pc_class_list));
    assertTrue("are subconstraint: " + pc_any + ", " + pc_class_arraylist, pc_any.subConstraintOf(pc_class_arraylist));

    // pc_class_string
    assertTrue("not subconstraint: " + pc_class_string + ", " + pc_any, !pc_class_string.subConstraintOf(pc_any));
    assertTrue("are subconstraint: " + pc_class_string + ", " + pc_class_string, pc_class_string.subConstraintOf(pc_class_string));
    assertTrue("not subconstraint: " + pc_class_string + ", " + pc_class_list, !pc_class_string.subConstraintOf(pc_class_list));
    assertTrue("not subconstraint: " + pc_class_string + ", " + pc_class_arraylist, !pc_class_string.subConstraintOf(pc_class_arraylist));

    // pc_class_list
    assertTrue("not subconstraint: " + pc_class_list + ", " + pc_any, !pc_class_list.subConstraintOf(pc_any));
    assertTrue("not subconstraint: " + pc_class_list + ", " + pc_class_string, !pc_class_list.subConstraintOf(pc_class_string));
    assertTrue("are subconstraint: " + pc_class_list + ", " + pc_class_list, pc_class_list.subConstraintOf(pc_class_list));
    assertTrue("are subconstraint: " + pc_class_list + ", " + pc_class_arraylist, pc_class_list.subConstraintOf(pc_class_arraylist));

    // pc_class_arraylist
    assertTrue("not subconstraint: " + pc_class_arraylist + ", " + pc_any, !pc_class_arraylist.subConstraintOf(pc_any));
    assertTrue("not subconstraint: " + pc_class_arraylist + ", " + pc_class_string, !pc_class_arraylist.subConstraintOf(pc_class_string));
    assertTrue("not subconstraint: " + pc_class_arraylist + ", " + pc_class_list, !pc_class_arraylist.subConstraintOf(pc_class_list));
    assertTrue("are subconstraint: " + pc_class_arraylist + ", " + pc_class_arraylist, pc_class_arraylist.subConstraintOf(pc_class_arraylist));
  }
  
  public void testExactProperties() {
    // self
    assertTrue("are subconstraints: " + pc_exact_fish + ", " + pc_exact_fish, pc_exact_fish.subConstraintOf(pc_exact_fish));
    assertTrue("are subconstraints: " + pc_exact_dog + ", " + pc_exact_dog, pc_exact_dog.subConstraintOf(pc_exact_dog));
    
    // not self
    assertTrue("not subconstraints: " + pc_exact_fish + ", " + pc_exact_dog, !pc_exact_fish.subConstraintOf(pc_exact_dog));
    assertTrue("not subconstraints: " + pc_exact_dog + ", " + pc_exact_fish, !pc_exact_dog.subConstraintOf(pc_exact_fish));
    
    // vs all
    assertTrue("not subconstraints: " + pc_exact_fish + ", " + pc_any, !pc_exact_fish.subConstraintOf(pc_any));
    assertTrue("are subconstraints: " + pc_any + ", " + pc_exact_fish, pc_any.subConstraintOf(pc_exact_fish));
    
    // vs pc_class_string
    assertTrue("not subconstraints: " + pc_exact_fish + ", " + pc_class_string, !pc_exact_fish.subConstraintOf(pc_class_string));
    assertTrue("are subconstraints: " + pc_class_string + ", " + pc_exact_fish, pc_class_string.subConstraintOf(pc_exact_fish));
    
    // vs pc_class_list
    assertTrue("not subconstraints: " + pc_exact_fish + ", " + pc_class_list, !pc_exact_fish.subConstraintOf(pc_class_list));
    assertTrue("not subconstraints: " + pc_class_list + ", " + pc_exact_fish, !pc_class_list.subConstraintOf(pc_exact_fish));
  }
  
  public void testEnumProperties() {
    // all
    assertTrue("are subconstraints: " + pc_any + ", " + pc_enum_animals, pc_any.subConstraintOf(pc_enum_animals));
    assertTrue("not subconstraints: " + pc_enum_animals + ", " + pc_any, !pc_enum_animals.subConstraintOf(pc_any));
    
    // compare them
    assertTrue("are subconstraints: " + pc_enum_animals + ", " + pc_enum_animals, pc_enum_animals.subConstraintOf(pc_enum_animals));
    assertTrue("are subconstraints: " + pc_enum_animals + ", " + pc_enum_birds, pc_enum_animals.subConstraintOf(pc_enum_birds));
    assertTrue("not subconstraints: " + pc_enum_animals + ", " + pc_enum_primary_colors, !pc_enum_animals.subConstraintOf(pc_enum_primary_colors));
    assertTrue("not subconstraints: " + pc_enum_birds + ", " + pc_enum_animals, !pc_enum_birds.subConstraintOf(pc_enum_animals));
    
    // to pc_class_string
    assertTrue("not subconstraints: " + pc_enum_animals + ", " + pc_class_string, !pc_enum_animals.subConstraintOf(pc_class_string));
    assertTrue("are subconstraints: " + pc_class_string + ", " + pc_enum_animals, pc_class_string.subConstraintOf(pc_enum_animals));
    
    // to pc_exact_dog
    assertTrue("are subconstraints: " + pc_enum_animals + ", " + pc_exact_dog, pc_enum_animals.subConstraintOf(pc_exact_dog));
    assertTrue("not subconstraints: " + pc_exact_dog + ", " + pc_enum_animals, !pc_exact_dog.subConstraintOf(pc_enum_animals));
    assertTrue("not subconstraints: " + pc_enum_primary_colors + ", " + pc_exact_dog, !pc_enum_primary_colors.subConstraintOf(pc_exact_dog));
  }
  
  public void testAllTypes() {
    // all
    assertTrue("are subconstraints: " + at_any + ", " + at_color, at_any.subTypeOf(at_color));
    assertTrue("not subconstraints: " + at_color + ", " + at_any, !at_color.subTypeOf(at_any));
  }
  
  public void testColorTypes() {
    assertTrue("are subtypes: " + at_colorsish + ", " + at_colorsish, at_colorsish.subTypeOf(at_colorsish));
    assertTrue("are subtypes: " + at_colorsish + ", " + at_colors, at_colorsish.subTypeOf(at_colors));
    assertTrue("are subtypes: " + at_colorsish + ", " + at_color, at_colorsish.subTypeOf(at_color));

    assertTrue("not subtypes: " + at_colors + ", " + at_colorsish, !at_colors.subTypeOf(at_colorsish));
    assertTrue("are subtypes: " + at_colors + ", " + at_colors, at_colors.subTypeOf(at_colors));
    assertTrue("are subtypes: " + at_colors + ", " + at_color, at_colors.subTypeOf(at_color));

    assertTrue("not subtypes: " + at_color + ", " + at_colorsish, !at_color.subTypeOf(at_colorsish));
    assertTrue("not subtypes: " + at_color + ", " + at_colors, !at_color.subTypeOf(at_colors));
    assertTrue("are subtypes: " + at_color + ", " + at_color, at_color.subTypeOf(at_color));
  }
  
  public void testColorName() {
    assertTrue("are subtypes: " + at_color + ", " + at_color, at_color.subTypeOf(at_color));
    assertTrue("are subtypes: " + at_color + ", " + at_color_name, at_color.subTypeOf(at_color_name));
    assertTrue("are subtypes: " + at_color + ", " + at_color, at_color_dog.subTypeOf(at_color_dog));

    assertTrue("not subtypes: " + at_color_name + ", " + at_color, !at_color_name.subTypeOf(at_color));
    assertTrue("are subtypes: " + at_color_name + ", " + at_color_name, at_color_name.subTypeOf(at_color_name));
    assertTrue("are subtypes: " + at_color_name + ", " + at_color_dog, at_color_name.subTypeOf(at_color_dog));

    assertTrue("not subtypes: " + at_color_dog + ", " + at_color, !at_color_dog.subTypeOf(at_color));
    assertTrue("not subtypes: " + at_color_dog + ", " + at_color_name, !at_color_dog.subTypeOf(at_color_name));
    assertTrue("are subtypes: " + at_color_dog + ", " + at_color_dog, at_color_dog.subTypeOf(at_color_dog));
  }
  
  public void testEnumerationOverlaps() {
      assertTrue(AnnotationTools.intersection(pc_enum_primary_colors, pc_enum_birds) == PropertyConstraint.NONE);
      assertTrue(AnnotationTools.intersection(pc_enum_primary_colors, pc_enum_warm_colors) != PropertyConstraint.NONE);
  }  
}
