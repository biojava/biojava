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
package org.biojava.utils.bytecode;

/**
 * Utility code for things you will frequently need.
 *
 * <p>
 * This class provides common constants representing access modifiers and
 * types. There are also some utility methods for munging data.
 * </p>
 *
 * @author Thomas Down
 * @author Matthew Pocock
 */
public class CodeUtils {
  public static final int ACC_PUBLIC       = 0x0001;
  public static final int ACC_PRIVATE      = 0x0002;
  public static final int ACC_PROTECTED    = 0x0004;
  public static final int ACC_STATIC       = 0x0008;
  public static final int ACC_FINAL        = 0x0010;
  public static final int ACC_SUPER        = 0x0020;
  public static final int ACC_SYNCHRONIZED = 0x0020;
  public static final int ACC_VOLATILE     = 0x0040;
  public static final int ACC_TRANSIENT    = 0x0080;
  public static final int ACC_NATIVE       = 0x0100;
  public static final int ACC_INTERFACE    = 0x0200;
  public static final int ACC_ABSTRACT     = 0x0400;
  public static final int ACC_STRICT       = 0x0800;
  
  public static final CodeClass TYPE_VOID;
  public static final CodeClass TYPE_INT;
  public static final CodeClass TYPE_FLOAT;
  public static final CodeClass TYPE_DOUBLE;
  public static final CodeClass TYPE_LONG;
  public static final CodeClass TYPE_BYTE;
  public static final CodeClass TYPE_SHORT;
  public static final CodeClass TYPE_CHAR;
  public static final CodeClass TYPE_BOOLEAN;
  public static final CodeClass TYPE_OBJECT;
  public static final CodeClass[] EMPTY_LIST;
  public static final CodeGenerator DO_NOTHING;
  
  static {
    TYPE_VOID = IntrospectedCodeClass.forClass(Void.TYPE);
    TYPE_BYTE = IntrospectedCodeClass.forClass(Byte.TYPE);
    TYPE_INT = IntrospectedCodeClass.forClass(Integer.TYPE);
    TYPE_FLOAT = IntrospectedCodeClass.forClass(Float.TYPE);
    TYPE_DOUBLE = IntrospectedCodeClass.forClass(Double.TYPE);
    TYPE_LONG = IntrospectedCodeClass.forClass(Long.TYPE);
    TYPE_SHORT = IntrospectedCodeClass.forClass(Short.TYPE);
    TYPE_CHAR = IntrospectedCodeClass.forClass(Character.TYPE);
    TYPE_BOOLEAN = IntrospectedCodeClass.forClass(Boolean.TYPE);
    TYPE_OBJECT = IntrospectedCodeClass.forClass(Object.class);
    
    EMPTY_LIST = new CodeClass[0];
    
    DO_NOTHING = new CodeGenerator() {
      public void writeCode(CodeContext cxt) { return; }
      public int stackDepth() { return 0; }
      public int stackDelta() { return 0; }
    };
  }
  
  /**
   * Format an array of classes as a comma-seperated list.
   *
   * <p>
   * The names of each class in classes will be seperated by a comma and a space
   * and will use CodeClass.getName() to produce strings for each one. Their
   * names will be present in the return value in the same order they are found
   * in the classes array
   * </p>
   *
   * @param classes  the array of classes to format
   * @return a String containing the list of class names
   */
  public static String classListToString(CodeClass[] classes) {
    StringBuffer sb = new StringBuffer();
    if(classes.length > 0) {
      sb.append(classes[0].getName());
    }
    
    for(int a = 1; a < classes.length; a++) {
      sb.append(", ");
      sb.append(classes[a].getName());
    }
    
    return sb.toString();
  }
  
  /**
   * Number of words needed for local variables of this type.
   *
   * <p>Longs and doubles require 2 words (64 bits), where as everything
   * else needs 1 word (32 bits). Void needs no words.
   * This just hides that knowledge.</p>
   *
   * @param cc the CodeClass to check word size for
   * @return number of words needed for this type
   */
  public static int wordsForType(CodeClass cc) {
    if(
      (cc == TYPE_DOUBLE) ||
      (cc == TYPE_LONG)
    ) {
      return 2;
    } else if(cc == TYPE_VOID) {
      return 0;
    } else {
      return 1;
    }
  }

  /**
   * Returns true if the class is a floating point number.
   *
   * <p>
   * Double and Float are floating point numbers. All other classes are not.
   * </p>
   *
   * @param cc  the class to check
   * @return  true if the class can be used to represent floating point numbers
   */
  public static boolean isFloatType(CodeClass cc) {
    return
      cc == TYPE_DOUBLE ||
      cc == TYPE_FLOAT;
  }

  /**
   * Returns true if the class is an integer number.
   *
   * <p>
   * All numeric types are integer (whole number) types, except for the floating
   * point types. All other classes are not integer types.
   * </p>
   *
   * @param cc  the class to check
   * @return  true if the class can be used to represent integer numbers
   */
  public static boolean isIntegerType(CodeClass cc) {
    return
      cc == TYPE_LONG ||
      cc == TYPE_BOOLEAN ||
      cc == TYPE_BYTE ||
      cc == TYPE_CHAR ||
      cc == TYPE_LONG ||
      cc == TYPE_SHORT;
  }
}
