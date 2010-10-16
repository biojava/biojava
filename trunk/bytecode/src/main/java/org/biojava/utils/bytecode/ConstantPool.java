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

import java.util.*;
import java.io.*;

/**
 * Build a Java class file constant pool.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 */

public class ConstantPool {
    private static final byte CONSTANT_Class = 7;
    private static final byte CONSTANT_Fieldref = 9;
    private static final byte CONSTANT_Methodref = 10;
    private static final byte CONSTANT_InterfaceMethodref = 11;
    private static final byte CONSTANT_String = 8;
    private static final byte CONSTANT_Integer = 3;
    private static final byte CONSTANT_Float = 4;
    private static final byte CONSTANT_Long = 5;
    private static final byte CONSTANT_Double = 6;
    private static final byte CONSTANT_NameAndType = 12;
    private static final byte CONSTANT_Utf8 = 1;

    private List constants;

    {
	constants = new ArrayList();
	constants.add(null); // Initial padder.
    }
   
    // Publically visible constant types

    public int resolveClass(CodeClass c) {
	return resolve(new CPTypedStringEntry(CONSTANT_Class, resolveUtf8(c.getJName())));
    }

    public int resolveField(CodeField f) {
      try {
        return resolve(new CPRefEntry(CONSTANT_Fieldref, resolveClass(f.getContainingClass()), resolveNameAndType(f.getName(), f.getType().getDescriptor())));
      } catch (NullPointerException npe) {
        throw new Error("Can't resolve filed " + f);
      }
    }

    public int resolveMethod(CodeMethod m) {
	return resolve(new CPRefEntry(CONSTANT_Methodref, resolveClass(m.getContainingClass()), resolveNameAndType(m.getName(), m.getDescriptor())));
    }

    public int resolveInterfaceMethod(CodeMethod m) {
	return resolve(new CPRefEntry(CONSTANT_InterfaceMethodref, resolveClass(m.getContainingClass()), resolveNameAndType(m.getName(), m.getDescriptor())));
    }

    public int resolveString(String s) {
	return resolve(new CPTypedStringEntry(CONSTANT_String, resolveUtf8(s)));
    }

    public int resolveInt(int i) {
	return resolve(new CPIntEntry(i));
    }

    public int resolveFloat(float f) {
	return resolve(new CPFloatEntry(f));
    }

    public int resolveLong(long l) {
	return resolve(new CPLongEntry(l));
    }

    public int resolveDouble(double d) {
	return resolve(new CPDoubleEntry(d));
    }

    // internal resolvers

    public int resolveUtf8(String s) {
	return resolve(new CPUtf8Entry(s));
    }

    public int resolveNameAndType(String name, String desc) {
	return resolve(new CPNameAndTypeEntry(resolveUtf8(name), resolveUtf8(desc)));
    }

    // The master resolver

    private int resolve(CPEntry e) {
	for (int i = 1; i < constants.size(); ++i) {
	    CPEntry e2 = (CPEntry) constants.get(i);
	    if (e2 != null && e.equals(e2))
		return i;
	}
	int i = constants.size();
	constants.add(e);
	for (int k = 1; k < e.needSlots(); ++k)
	    constants.add(null);
	return i;
    }

    // Output again

    public int constantPoolSize() {
	return constants.size();
    }

  public void writeConstantPool(DataOutput d) throws IOException {
//    int count = 1;
    for (Iterator i = constants.iterator(); i.hasNext(); ) {
	    CPEntry e = (CPEntry) i.next();
	    if (e != null) {
//        System.out.println("Writing constant " + count + " " + e);
//        count += e.needSlots();
        e.write(d);
      }
    }
  }
  
    // Types for storing the cpool

    private static interface CPEntry {
	public void write(DataOutput d) throws IOException;
	public int needSlots();
    }

    private static class CPTypedStringEntry implements CPEntry {
	byte tag;
	int name;

        CPTypedStringEntry(byte tag, int name) {
	    this.tag = tag;
	    this.name = name;
	}

	public boolean equals(Object o) {
	    if (! (o instanceof CPTypedStringEntry))
		return false;

	    CPTypedStringEntry cte = (CPTypedStringEntry) o;
	    return (cte.name == name && cte.tag == tag);
	}

	public void write(DataOutput d) throws IOException {
	    d.writeByte(tag);
	    d.writeShort(name);
	}

	public int needSlots() {
	    return 1;
	}
  
  public String toString() {
    return "CPTypedStringEntry tag=" + tag + " name=" + name;
  }
    }

    private static class CPRefEntry implements CPEntry {
	byte tag;
	int clazz;
	int name;

        CPRefEntry(byte tag, int clazz, int name) {
	    this.tag = tag;
	    this.clazz = clazz;
	    this.name = name;
	}

	public boolean equals(Object o) {
	    if (! (o instanceof CPRefEntry))
		return false;

	    CPRefEntry cte = (CPRefEntry) o;
	    return (cte.clazz == clazz && cte.name == name && cte.tag == tag);
	}

	public void write(DataOutput d) throws IOException {
	    d.writeByte(tag);
	    d.writeShort(clazz);
	    d.writeShort(name);
	}

	public int needSlots() {
	    return 1;
	}
  
  public String toString() {
    return "CPRefEntry tag=" + tag + " class=" + clazz + " name=" + name;
  }
    }

    private static class CPIntEntry implements CPEntry {
	int val;

        CPIntEntry(int val) {
	    this.val = val;
	}

	public boolean equals(Object o) {
	    if (! (o instanceof CPIntEntry))
		return false;

	    return (((CPIntEntry) o).val == val);
	}

	public void write(DataOutput d) throws IOException {
	    d.writeByte(CONSTANT_Integer);
	    d.writeInt(val);
	}

	public int needSlots() {
	    return 1;
	}
  
  public String toString() {
    return "CPIntEntry val=" + val;
  }
    }

    private static class CPLongEntry implements CPEntry {
	long val;

        CPLongEntry(long val) {
	    this.val = val;
	}

	public boolean equals(Object o) {
	    if (! (o instanceof CPLongEntry))
		return false;

	    return (((CPLongEntry) o).val == val);
	}

	public void write(DataOutput d) throws IOException {
	    d.writeByte(CONSTANT_Long);
	    d.writeLong(val);
	}

	public int needSlots() {
	    return 2;
	}
  
  public String toString() {
    return "CPLongEntry val=" + val;
  }
    }

    private static class CPFloatEntry implements CPEntry {
	float val;

        CPFloatEntry(float val) {
	    this.val = val;
	}

	public boolean equals(Object o) {
	    if (! (o instanceof CPFloatEntry))
		return false;

	    return
        (((CPFloatEntry) o).val == val) ||
        (Float.isNaN(((CPFloatEntry) o).val) && Float.isNaN(val));
	}

	public void write(DataOutput d) throws IOException {
	    d.writeByte(CONSTANT_Float);
	    d.writeFloat(val);
	}

	public int needSlots() {
	    return 1;
	}
  
  public String toString() {
    return "CPFloatEntry val=" + val;
  }
    }

    private static class CPDoubleEntry implements CPEntry {
	double val;

        CPDoubleEntry(double val) {
	    this.val = val;
	}

	public boolean equals(Object o) {
	    if (! (o instanceof CPDoubleEntry))
		return false;

	    return
        (((CPDoubleEntry) o).val == val) ||
        (Double.isNaN(((CPDoubleEntry) o).val) && Double.isNaN(val));
	}

	public void write(DataOutput d) throws IOException {
	    d.writeByte(CONSTANT_Double);
	    d.writeDouble(val);
	}

	public int needSlots() {
	    return 2;
	}
  
  public String toString() {
    return "CPDoubleEntry val=" + val;
  }
    }

    private static class CPUtf8Entry implements CPEntry {
	String val;

        CPUtf8Entry(String val) {
	    this.val = val;
	}

	public boolean equals(Object o) {
	    if (! (o instanceof CPUtf8Entry))
		return false;

	    return (((CPUtf8Entry) o).val.equals(val));
	}

	public void write(DataOutput d) throws IOException {
	    d.writeByte(CONSTANT_Utf8);
	    d.writeUTF(val);
	}

	public int needSlots() {
	    return 1;
	}
  
  public String toString() {
    return "CPUtf8Entry val=" + val;
  }
    }

    private static class CPNameAndTypeEntry implements CPEntry {
	int name;
	int desc;

        CPNameAndTypeEntry(int name, int desc) {
	    this.name = name;
	    this.desc = desc;
	}

	public boolean equals(Object o) {
	    if (! (o instanceof CPNameAndTypeEntry))
		return false;

	    CPNameAndTypeEntry cpnte = (CPNameAndTypeEntry) o;
	    return (cpnte.desc == desc && cpnte.name == name);
	}

	public void write(DataOutput d) throws IOException {
	    d.writeByte(CONSTANT_NameAndType);
	    d.writeShort(name);
	    d.writeShort(desc);
	}

	public int needSlots() {
	    return 1;
	}
  
  public String toString() {
    return "CPNameAndTypeEntry name=" + name + " desc=" + desc;
  }
    }
}
