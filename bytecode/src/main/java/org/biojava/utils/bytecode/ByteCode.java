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
 * Factory for objects which encapsulate individual Java bytecode instructions.
 * Most methods in this class are auto-generated.
 *
 * <p>There are two classes of methods. The first ones are for creating objects
 * that directly represent opcodes. These effectively wrap one of the opcode
 * constants. The others do something more clever. For example, make_if emits
 * something that is equivalent to a normal Java if statement.</p>
 *
 * <p>Generic types are supported using the factory methods that take
 * ParametricType arguments.
 *
 * @author Thomas Down
 * @author Matthew Pocock
 */

public class ByteCode {
    //
    // iconst and friends
    //

    public static Instruction make_iconst(int i) {
      if (i >= -1 && i <= 5) {
        return new NoOperandsInstruction((byte) (op_iconst_0 + i), 1);
      } else if (i >= Byte.MIN_VALUE && i <= Byte.MAX_VALUE) {
        return new ByteInstruction(op_bipush, (byte) i, 1);
      } else if (i >= Short.MIN_VALUE && i <= Short.MAX_VALUE) {
        return new ShortInstruction(op_sipush, i, 1);
      }
      
      return new IntConstantInstruction(i);
    }

    //
    // String constants
    //

    public static Instruction make_sconst(String s) {
	return new StringConstantInstruction(s);
    }

    //
    // Double constants
    //

    public static Instruction make_dconst(double d) {
      if (d == 0.0) {
        return new NoOperandsInstruction(op_dconst_0, 1);
      } else if (d == 1.0) {
        return new NoOperandsInstruction(op_dconst_1, 1);
      } else {
        return new DoubleInstruction(d);
      }
    }

    //
    // Long constants
    //

    public static Instruction make_lconst(long l) {
      if (l == 0L) {
        return new NoOperandsInstruction(op_lconst_0, 1);
      } else if (l == 1L) {
        return new NoOperandsInstruction(op_lconst_1, 1);
      } else {
        return new LongConstantInstruction(l);
      }
    }

    //
    // Float constants
    //

    public static Instruction make_fconst(float f) {
      if (f == 0.0F) {
        return new NoOperandsInstruction(op_fconst_0, 1);
      } else if (f == 1.0F) {
        return new NoOperandsInstruction(op_fconst_1, 1);
      } else if (f == 2.0F) {
        return new NoOperandsInstruction(op_fconst_2, 1);
      } else {
        return new FloatConstantInstruction(f);
      }
    }

    //
    // Generate instruction objects for invokes
    //

    public static Instruction make_invokevirtual(CodeMethod cm) {
	return new MethodInstruction(op_invokevirtual, cm);
    }

    public static Instruction make_invokespecial(CodeMethod cm) {
	return new MethodInstruction(op_invokespecial, cm);
    }

    public static Instruction make_invokestatic(CodeMethod cm) {
	return new MethodInstruction(op_invokestatic, cm);
    }

    public static Instruction make_invokeinterface(CodeMethod cm) {
	return new MethodInstruction(op_invokeinterface, cm);
    }

    //
    // Generate instruction objects for fields
    //

    public static Instruction make_getfield(CodeField cf) {
	return new FieldInstruction(op_getfield, cf);
    }

    public static Instruction make_putfield(CodeField cf) {
	return new FieldInstruction(op_putfield, cf);
    }

    public static Instruction make_getstatic(CodeField cf) {
	return new FieldInstruction(op_getstatic, cf);
    }

    public static Instruction make_putstatic(CodeField cf) {
	return new FieldInstruction(op_putstatic, cf);
    }

    //
    // Generate instruction objects for loads
    //

    public static Instruction make_iload(LocalVariable lv) 
        throws CodeException
    {
	CodeClass cc = lv.getType();
	if (cc != CodeUtils.TYPE_INT && 
	    cc != CodeUtils.TYPE_SHORT && 
	    cc != CodeUtils.TYPE_CHAR &&
	    cc != CodeUtils.TYPE_BYTE &&
	    cc != CodeUtils.TYPE_BOOLEAN)
	{
	    throw new CodeException(cc.getName() + " is not a VM `i' type");
	}

	return new LocalVariableInstruction(op_iload, op_iload_0, lv);
    }

    public static Instruction make_lload(LocalVariable lv) 
        throws CodeException
    {
	if (lv.getType() != CodeUtils.TYPE_LONG) {
	    throw new CodeException(lv.getType().getName() + " is not a long");
	}

	return new LocalVariableInstruction(op_lload, op_lload_0, lv);
    }

    public static Instruction make_fload(LocalVariable lv) 
        throws CodeException
    {
	if (lv.getType() != CodeUtils.TYPE_FLOAT) {
	    throw new CodeException(lv.getType().getName() + " is not a float");
	}

	return new LocalVariableInstruction(op_fload, op_fload_0, lv);
    }

    public static Instruction make_dload(LocalVariable lv) 
        throws CodeException
    {
	if (lv.getType() != CodeUtils.TYPE_DOUBLE) {
	    throw new CodeException(lv.getType().getName() + " is not a double");
	}

	return new LocalVariableInstruction(op_dload, op_dload_0, lv);
    }

    public static Instruction make_aload(LocalVariable lv) 
        throws CodeException
    {
	if (lv.getType().isPrimitive()) {
	    throw new CodeException(lv.getType().getName() + " is a primitive type");
	}

	return new LocalVariableInstruction(op_aload, op_aload_0, lv);
    }

    //
    // Generate instruction objects for stores
    //

    public static Instruction make_istore(LocalVariable lv) 
        throws CodeException
    {
	CodeClass cc = lv.getType();
	if (cc != CodeUtils.TYPE_INT && 
	    cc != CodeUtils.TYPE_SHORT && 
	    cc != CodeUtils.TYPE_CHAR &&
	    cc != CodeUtils.TYPE_BYTE &&
	    cc != CodeUtils.TYPE_BOOLEAN)
	{
	    throw new CodeException(cc.getName() + " is not a VM `i' type");
	}


	return new LocalVariableInstruction(op_istore, op_istore_0, lv);
    }

    public static Instruction make_lstore(LocalVariable lv)
        throws CodeException
    {
	if (lv.getType() != CodeUtils.TYPE_LONG) {
	    throw new CodeException(lv.getType().getName() + " is not a long");
	}

	return new LocalVariableInstruction(op_lstore, op_lstore_0, lv);
    }

    public static Instruction make_fstore(LocalVariable lv) 
        throws CodeException
    {
	if (lv.getType() != CodeUtils.TYPE_FLOAT) {
	    throw new CodeException(lv.getType().getName() + " is not a float");
	}

	return new LocalVariableInstruction(op_fstore, op_fstore_0, lv);
    }

    public static Instruction make_dstore(LocalVariable lv) 
        throws CodeException
    {
	if (lv.getType() != CodeUtils.TYPE_DOUBLE) {
	    throw new CodeException(lv.getType().getName() + " is not a double");
	}

	return new LocalVariableInstruction(op_dstore, op_dstore_0, lv);
    }

    public static Instruction make_astore(LocalVariable lv) 
        throws CodeException
    {
	if (lv.getType().isPrimitive()) {
	    throw new CodeException(lv.getType().getName() + " is a primitive type");
	}

	return new LocalVariableInstruction(op_astore, op_astore_0, lv);
    }

    //
    // Generate the Instruction objects for all this `if's
    //

    /**
     * Return the Instruction object for the ifeq instruction.
     */

    public static Instruction make_ifeq(Label lab) {
        return make_if(op_ifeq, lab);
    }

    /**
     * Return the Instruction object for the ifne instruction.
     */

    public static Instruction make_ifne(Label lab) {
        return make_if(op_ifne, lab);
    }

    /**
     * Return the Instruction object for the iflt instruction.
     *
     */

    public static Instruction make_iflt(Label lab) {
        return new LabelInstruction(op_iflt, lab, -1);
    }

    /**
     * Return the Instruction object for the ifge instruction.
     *
     */

    public static Instruction make_ifge(Label lab) {
        return new LabelInstruction(op_ifge, lab, -1);
    }

    /**
     * Return the Instruction object for the ifgt instruction.
     *
     */

    public static Instruction make_ifgt(Label lab) {
        return new LabelInstruction(op_ifgt, lab, -1);
    }

    /**
     * Return the Instruction object for the ifle instruction.
     */

    public static Instruction make_ifle(Label lab) {
        return new LabelInstruction(op_ifle, lab, -1);
    }

    /**
     * Return the Instruction object for the if_icmpeq instruction.
     */

    public static Instruction make_if_icmpeq(Label l) {
        return new LabelInstruction(op_if_icmpeq, l, -2);
    }

    /**
     * Return the Instruction object for the if_icmpne instruction.
     */

    public static Instruction make_if_icmpne(Label l) {
        return new LabelInstruction(op_if_icmpne, l, -2);
    }

    /**
     * Return the Instruction object for the if_icmplt instruction.
     */

    public static Instruction make_if_icmplt(Label l) {
        return new LabelInstruction(op_if_icmplt, l, -2);
    }

    /**
     * Return the Instruction object for the if_icmpge instruction.
     */

    public static Instruction make_if_icmpge(Label l) {
        return new LabelInstruction(op_if_icmpge, l, -2);
    }

    /**
     * Return the Instruction object for the if_icmpgt instruction.
     */

    public static Instruction make_if_icmpgt(Label l) {
        return new LabelInstruction(op_if_icmpgt, l, -2);
    }

    /**
     * Return the Instruction object for the if_icmple instruction.
     */

    public static Instruction make_if_icmple(Label l) {
        return new LabelInstruction(op_if_icmple, l, -2);
    }

    /**
     * Return the Instruction object for the if_acmpeq instruction.
     */

    public static Instruction make_if_acmpeq(Label l) {
        return new LabelInstruction(op_if_acmpeq, l, -2);
    }

    /**
     * Return the Instruction object for the if_acmpne instruction.
     */

    public static Instruction make_if_acmpne(Label l) {
        return new LabelInstruction(op_if_acmpne, l, -2);
    }

    /**
     * Return the Instruction object for the ifnull instruction.
     */

    public static Instruction make_ifnull(Label l) {
        return new LabelInstruction(op_ifnull, l, -1);
    }

    /**
     * Return the Instruction object for the ifnonnull instruction.
     */

    public static Instruction make_ifnonnull(Label l) {
        return new LabelInstruction(op_ifnonnull, l, -1);
    }

    //
    // The other label instructions
    //

    public static Instruction make_goto(Label l) {
        return new LabelInstruction(op_goto, l, 0);
    }

    public static Instruction make_jsr(Label l) {
        return new LabelInstruction(op_jsr, l, 1);
    }

    //
    // Generate the Instruction objects for all the no-operands instructions
    //

    /**
     * Return the Instruction object for the nop instruction.
     */

    public static Instruction make_nop() {
        return new NoOperandsInstruction(op_nop, 0);
    }

    /**
     * Return the Instruction object for the aconst_null instruction.
     */

    public static Instruction make_aconst_null() {
        return new NoOperandsInstruction(op_aconst_null, 1);
    }

    /**
     * Return the Instruction object for the iaload instruction.
     */

    public static Instruction make_iaload() {
        return new NoOperandsInstruction(op_iaload, -1);
    }

    /**
     * Return the Instruction object for the laload instruction.
     */

    public static Instruction make_laload() {
        return new NoOperandsInstruction(op_laload, -1);
    }

    /**
     * Return the Instruction object for the faload instruction.
     */

    public static Instruction make_faload() {
        return new NoOperandsInstruction(op_faload, -1);
    }

    /**
     * Return the Instruction object for the daload instruction.
     */

    public static Instruction make_daload() {
        return new NoOperandsInstruction(op_daload, -1);
    }

    /**
     * Return the Instruction object for the aaload instruction.
     */

    public static Instruction make_aaload() {
        return new NoOperandsInstruction(op_aaload, -1);
    }

    /**
     * Return the Instruction object for the baload instruction.
     */

    public static Instruction make_baload() {
        return new NoOperandsInstruction(op_baload, -1);
    }

    /**
     * Return the Instruction object for the caload instruction.
     */

    public static Instruction make_caload() {
        return new NoOperandsInstruction(op_caload, -1);
    }

    /**
     * Return the Instruction object for the saload instruction.
     */

    public static Instruction make_saload() {
        return new NoOperandsInstruction(op_saload, -1);
    }

    /**
     * Return the Instruction object for the iastore instruction.
     */

    public static Instruction make_iastore() {
        return new NoOperandsInstruction(op_iastore, -3);
    }

    /**
     * Return the Instruction object for the lastore instruction.
     */

    public static Instruction make_lastore() {
        return new NoOperandsInstruction(op_lastore, -3);
    }

    /**
     * Return the Instruction object for the fastore instruction.
     */

    public static Instruction make_fastore() {
        return new NoOperandsInstruction(op_fastore, -3);
    }

    /**
     * Return the Instruction object for the dastore instruction.
     */

    public static Instruction make_dastore() {
        return new NoOperandsInstruction(op_dastore, -3);
    }

    /**
     * Return the Instruction object for the aastore instruction.
     */

    public static Instruction make_aastore() {
        return new NoOperandsInstruction(op_aastore, -3);
    }

    /**
     * Return the Instruction object for the bastore instruction.
     */

    public static Instruction make_bastore() {
        return new NoOperandsInstruction(op_bastore, -3);
    }

    /**
     * Return the Instruction object for the castore instruction.
     */

    public static Instruction make_castore() {
        return new NoOperandsInstruction(op_castore, -3);
    }

    /**
     * Return the Instruction object for the sastore instruction.
     */

    public static Instruction make_sastore() {
        return new NoOperandsInstruction(op_sastore, -3);
    }

    /**
     * Return the Instruction object for the pop instruction.
     */

    public static Instruction make_pop() {
        return new NoOperandsInstruction(op_pop, -1);
    }

    /**
     * Return the Instruction object for the pop2 instruction.
     */

    public static Instruction make_pop2() {
        return new NoOperandsInstruction(op_pop2, -2);
    }

    /**
     * Return the Instruction object for the dup instruction.
     */

    public static Instruction make_dup() {
        return new NoOperandsInstruction(op_dup, 1);
    }

    /**
     * Return the Instruction object for the dup_x1 instruction.
     */

    public static Instruction make_dup_x1() {
        return new NoOperandsInstruction(op_dup_x1, 1);
    }

    /**
     * Return the Instruction object for the dup_x2 instruction.
     */

    public static Instruction make_dup_x2() {
        return new NoOperandsInstruction(op_dup_x2, 1);
    }

    /**
     * Return the Instruction object for the dup2 instruction.
     */

    public static Instruction make_dup2() {
      return new NoOperandsInstruction(op_dup2, 2);
    }

    /**
     * Return the Instruction object for the dup2_x1 instruction.
     */

    public static Instruction make_dup2_x1() {
      return new NoOperandsInstruction(op_dup2_x1, 2);
    }

    /**
     * Return the Instruction object for the dup2_x2 instruction.
     */

    public static Instruction make_dup2_x2() {
        return new NoOperandsInstruction(op_dup2_x2, 2);
    }

    /**
     * Return the Instruction object for the swap instruction.
     */

    public static Instruction make_swap() {
        return new NoOperandsInstruction(op_swap, 0);
    }

    /**
     * Return the Instruction object for the iadd instruction.
     */

    public static Instruction make_iadd() {
        return new NoOperandsInstruction(op_iadd, -1);
    }

    /**
     * Return the Instruction object for the ladd instruction.
     */

    public static Instruction make_ladd() {
        return new NoOperandsInstruction(op_ladd, -1);
    }

    /**
     * Return the Instruction object for the fadd instruction.
     */

    public static Instruction make_fadd() {
        return new NoOperandsInstruction(op_fadd, -1);
    }

    /**
     * Return the Instruction object for the dadd instruction.
     */

    public static Instruction make_dadd() {
        return new NoOperandsInstruction(op_dadd, -1);
    }

    /**
     * Return the Instruction object for the isub instruction.
     */

    public static Instruction make_isub() {
        return new NoOperandsInstruction(op_isub, -1);
    }

    /**
     * Return the Instruction object for the lsub instruction.
     */

    public static Instruction make_lsub() {
        return new NoOperandsInstruction(op_lsub, -1);
    }

    /**
     * Return the Instruction object for the fsub instruction.
     */

    public static Instruction make_fsub() {
        return new NoOperandsInstruction(op_fsub, -1);
    }

    /**
     * Return the Instruction object for the dsub instruction.
     */

    public static Instruction make_dsub() {
        return new NoOperandsInstruction(op_dsub, -1);
    }

    /**
     * Return the Instruction object for the imul instruction.
     */

    public static Instruction make_imul() {
        return new NoOperandsInstruction(op_imul, -1);
    }

    /**
     * Return the Instruction object for the lmul instruction.
     */

    public static Instruction make_lmul() {
        return new NoOperandsInstruction(op_lmul, -1);
    }

    /**
     * Return the Instruction object for the fmul instruction.
     */

    public static Instruction make_fmul() {
        return new NoOperandsInstruction(op_fmul, -1);
    }

    /**
     * Return the Instruction object for the dmul instruction.
     */

    public static Instruction make_dmul() {
        return new NoOperandsInstruction(op_dmul, -1);
    }

    /**
     * Return the Instruction object for the idiv instruction.
     */

    public static Instruction make_idiv() {
        return new NoOperandsInstruction(op_idiv, -1);
    }

    /**
     * Return the Instruction object for the ldiv instruction.
     */

    public static Instruction make_ldiv() {
        return new NoOperandsInstruction(op_ldiv, -1);
    }

    /**
     * Return the Instruction object for the fdiv instruction.
     */

    public static Instruction make_fdiv() {
        return new NoOperandsInstruction(op_fdiv, -1);
    }

    /**
     * Return the Instruction object for the ddiv instruction.
     */

    public static Instruction make_ddiv() {
        return new NoOperandsInstruction(op_ddiv, -1);
    }

    /**
     * Return the Instruction object for the irem instruction.
     */

    public static Instruction make_irem() {
        return new NoOperandsInstruction(op_irem, -1);
    }

    /**
     * Return the Instruction object for the lrem instruction.
     */

    public static Instruction make_lrem() {
        return new NoOperandsInstruction(op_lrem, -1);
    }

    /**
     * Return the Instruction object for the frem instruction.
     */

    public static Instruction make_frem() {
        return new NoOperandsInstruction(op_frem, -1);
    }

    /**
     * Return the Instruction object for the drem instruction.
     */

    public static Instruction make_drem() {
        return new NoOperandsInstruction(op_drem, -1);
    }

    /**
     * Return the Instruction object for the ineg instruction.
     */

    public static Instruction make_ineg() {
        return new NoOperandsInstruction(op_ineg, 0);
    }

    /**
     * Return the Instruction object for the lneg instruction.
     */

    public static Instruction make_lneg() {
        return new NoOperandsInstruction(op_lneg, 0);
    }

    /**
     * Return the Instruction object for the fneg instruction.
     */

    public static Instruction make_fneg() {
        return new NoOperandsInstruction(op_fneg, 0);
    }

    /**
     * Return the Instruction object for the dneg instruction.
     */

    public static Instruction make_dneg() {
        return new NoOperandsInstruction(op_dneg, 0);
    }

    /**
     * Return the Instruction object for the ishl instruction.
     */

    public static Instruction make_ishl() {
        return new NoOperandsInstruction(op_ishl, -1);
    }

    /**
     * Return the Instruction object for the lshl instruction.
     */

    public static Instruction make_lshl() {
        return new NoOperandsInstruction(op_lshl, -1);
    }

    /**
     * Return the Instruction object for the ishr instruction.
     */

    public static Instruction make_ishr() {
        return new NoOperandsInstruction(op_ishr, -1);
    }

    /**
     * Return the Instruction object for the lshr instruction.
     */

    public static Instruction make_lshr() {
        return new NoOperandsInstruction(op_lshr, -1);
    }

    /**
     * Return the Instruction object for the iushr instruction.
     */

    public static Instruction make_iushr() {
        return new NoOperandsInstruction(op_iushr, -1);
    }

    /**
     * Return the Instruction object for the lushr instruction.
     */

    public static Instruction make_lushr() {
        return new NoOperandsInstruction(op_lushr, -1);
    }

    /**
     * Return the Instruction object for the iand instruction.
     */

    public static Instruction make_iand() {
        return new NoOperandsInstruction(op_iand, -1);
    }

    /**
     * Return the Instruction object for the land instruction.
     */

    public static Instruction make_land() {
        return new NoOperandsInstruction(op_land, -1);
    }

    /**
     * Return the Instruction object for the ior instruction.
     */

    public static Instruction make_ior() {
        return new NoOperandsInstruction(op_ior, -1);
    }

    /**
     * Return the Instruction object for the lor instruction.
     */

    public static Instruction make_lor() {
        return new NoOperandsInstruction(op_lor, -1);
    }

    /**
     * Return the Instruction object for the ixor instruction.
     */

    public static Instruction make_ixor() {
        return new NoOperandsInstruction(op_ixor, -1);
    }

    /**
     * Return the Instruction object for the lxor instruction.
     */

    public static Instruction make_lxor() {
        return new NoOperandsInstruction(op_lxor, -1);
    }

    /**
     * Return the Instruction object for the i2l instruction.
     */

    public static Instruction make_i2l() {
        return new NoOperandsInstruction(op_i2l, 0);
    }

    /**
     * Return the Instruction object for the i2f instruction.
     */

    public static Instruction make_i2f() {
        return new NoOperandsInstruction(op_i2f, 0);
    }

    /**
     * Return the Instruction object for the i2d instruction.
     */

    public static Instruction make_i2d() {
        return new NoOperandsInstruction(op_i2d, 0);
    }

    /**
     * Return the Instruction object for the l2i instruction.
     */

    public static Instruction make_l2i() {
        return new NoOperandsInstruction(op_l2i, 0);
    }

    /**
     * Return the Instruction object for the l2f instruction.
     */

    public static Instruction make_l2f() {
        return new NoOperandsInstruction(op_l2f, 0);
    }

    /**
     * Return the Instruction object for the l2d instruction.
     */

    public static Instruction make_l2d() {
        return new NoOperandsInstruction(op_l2d, 0);
    }

    /**
     * Return the Instruction object for the f2i instruction.
     */

    public static Instruction make_f2i() {
        return new NoOperandsInstruction(op_f2i, 0);
    }

    /**
     * Return the Instruction object for the f2l instruction.
     */

    public static Instruction make_f2l() {
        return new NoOperandsInstruction(op_f2l, 0);
    }

    /**
     * Return the Instruction object for the f2d instruction.
     */

    public static Instruction make_f2d() {
        return new NoOperandsInstruction(op_f2d, 0);
    }

    /**
     * Return the Instruction object for the d2i instruction.
     */

    public static Instruction make_d2i() {
        return new NoOperandsInstruction(op_d2i, 0);
    }

    /**
     * Return the Instruction object for the d2l instruction.
     */

    public static Instruction make_d2l() {
        return new NoOperandsInstruction(op_d2l, 0);
    }

    /**
     * Return the Instruction object for the d2f instruction.
     */

    public static Instruction make_d2f() {
        return new NoOperandsInstruction(op_d2f, 0);
    }

    /**
     * Return the Instruction object for the i2b instruction.
     */

    public static Instruction make_i2b() {
        return new NoOperandsInstruction(op_i2b, 0);
    }

    /**
     * Return the Instruction object for the i2c instruction.
     */

    public static Instruction make_i2c() {
        return new NoOperandsInstruction(op_i2c, 0);
    }

    /**
     * Return the Instruction object for the i2s instruction.
     */

    public static Instruction make_i2s() {
        return new NoOperandsInstruction(op_i2s, 0);
    }

    /**
     * Return the Instruction object for the lcmp instruction.
     */

    public static Instruction make_lcmp() {
        return new NoOperandsInstruction(op_lcmp, -1);
    }

    /**
     * Return the Instruction object for the fcmpl instruction.
     */

    public static Instruction make_fcmpl() {
        return new NoOperandsInstruction(op_fcmpl, -1);
    }

    /**
     * Return the Instruction object for the fcmpg instruction.
     */

    public static Instruction make_fcmpg() {
        return new NoOperandsInstruction(op_fcmpg, -1);
    }

    /**
     * Return the Instruction object for the dcmpl instruction.
     */

    public static Instruction make_dcmpl() {
        return new NoOperandsInstruction(op_dcmpl, -1);
    }

    /**
     * Return the Instruction object for the dcmpg instruction.
     */

    public static Instruction make_dcmpg() {
        return new NoOperandsInstruction(op_dcmpg, -1);
    }

    /**
     * Return the Instruction object for the ireturn instruction.
     */

    public static Instruction make_ireturn() {
        return new NoOperandsInstruction(op_ireturn, -1);
    }

    /**
     * Return the Instruction object for the lreturn instruction.
     */

    public static Instruction make_lreturn() {
        return new NoOperandsInstruction(op_lreturn, -1);
    }

    /**
     * Return the Instruction object for the freturn instruction.
     */

    public static Instruction make_freturn() {
        return new NoOperandsInstruction(op_freturn, -1);
    }

    /**
     * Return the Instruction object for the dreturn instruction.
     */

    public static Instruction make_dreturn() {
        return new NoOperandsInstruction(op_dreturn, -1);
    }

    /**
     * Return the Instruction object for the areturn instruction.
     */

    public static Instruction make_areturn() {
        return new NoOperandsInstruction(op_areturn, -1);
    }

    /**
     * Return the Instruction object for the return instruction.
     */

    public static Instruction make_return() {
        return new NoOperandsInstruction(op_return, 0);
    }
    
    /**
     * Return the Instruction object for the arraylength instruction.
     */

    public static Instruction make_arraylength() {
        return new NoOperandsInstruction(op_arraylength, 0);
    }

    /**
     * Return the Instruction object for the athrow instruction.
     */

    public static Instruction make_athrow() {
        return new NoOperandsInstruction(op_athrow, 0);
    }

    /**
     * Return the Instruction object for the monitorenter instruction.
     *
     */

    public static Instruction make_monitorenter() {
        return new NoOperandsInstruction(op_monitorenter, -1);
    }

    /**
     * Return the Instruction object for the monitorexit instruction.
     */

    public static Instruction make_monitorexit() {
        return new NoOperandsInstruction(op_monitorexit, -1);
    }

    /**
     * Return the Instruction object for the wide instruction.
     */

    public static Instruction make_wide() {
        return new NoOperandsInstruction(op_wide, 0);
    }

    /**
     * Return the Instruction object for the breakpoint instruction.
     */

    public static Instruction make_breakpoint() {
        return new NoOperandsInstruction(op_breakpoint, 0);
    }

    /**
     * Return the Instruction object for the invokevirtual_quick_w instruction.
     *
     * Apparently now undocumented. Presumably deprecated.
     */

    public static Instruction make_invokevirtual_quick_w() {
        // fixme - I have no idea what the stack-depth change is here
        // should we even be using this class?
        return new NoOperandsInstruction(op_invokevirtual_quick_w, 0);
    }

    /**
     * Return the Instruction object for the impdep1 instruction.
     */

    public static Instruction make_impdep1() {
        return new NoOperandsInstruction(op_impdep1, 0);
    }

    /**
     * Return the Instruction object for the impdep2 instruction.
     */

    public static Instruction make_impdep2() {
        return new NoOperandsInstruction(op_impdep2, 0);
    }

  public static Instruction make_checkcast(CodeClass clazz) {
    return new ClassInstruction(op_checkcast, clazz, 0);
  }

  public static Instruction make_instanceof(CodeClass clazz) {
    return new ClassInstruction(op_instanceof, clazz, 0);
  }
    // new instruction added by hand - mrp
    
    public static Instruction make_new(CodeClass clazz) {
      return new ClassInstruction(op_new, clazz, 1);
    }

    public static Instruction make_newarray(CodeClass clazz) 
        throws CodeException
    {
      if (clazz.isPrimitive()) {
        int type = -1;
        if (clazz == CodeUtils.TYPE_BOOLEAN) {
          type = 4;
        } else if (clazz == CodeUtils.TYPE_CHAR) {
          type = 5;
        } else if (clazz == CodeUtils.TYPE_FLOAT) {
          type = 6;
        } else if (clazz == CodeUtils.TYPE_DOUBLE) {
          type = 7;
        } else if (clazz == CodeUtils.TYPE_BYTE) {
          type = 8;
        } else if (clazz == CodeUtils.TYPE_SHORT) {
          type = 9;
        } else if (clazz == CodeUtils.TYPE_INT) {
          type = 10;
        } else if (clazz == CodeUtils.TYPE_LONG) {
          type = 11;
        } 
        if (type < 0) { 
          throw new CodeException("Invalid type " + clazz.getName());
        }
       
        return new ByteInstruction(op_newarray, (byte) type, 0);
      } else {
        return new ClassInstruction(op_anewarray, clazz, 0);
      }
    }

    /**
     * A convenient one-stop method to get a return statement suitable for a
     * method.
     *
     * @param method  the CodeMethod to return from
     * @return a return Instruction that fits this method's return type
     */
    public static Instruction make_return(CodeMethod method) {
      return make_return(method.getReturnType());
    }
    
    /**
     * Creates the return Instruction suitable for a given class or type.
     *
     * @param clazz  the CodeClass representing the class or type to return
     * @return  the apropreate return instruction
     */
    public static Instruction make_return(CodeClass clazz) {
      if (false) {
      } else if(CodeUtils.TYPE_VOID.equals(clazz)) {
        return make_return();
      } else if(
        CodeUtils.TYPE_BYTE.equals(clazz) ||
        CodeUtils.TYPE_SHORT.equals(clazz) ||
        CodeUtils.TYPE_CHAR.equals(clazz) ||
        CodeUtils.TYPE_BOOLEAN.equals(clazz) ||
        CodeUtils.TYPE_INT.equals(clazz)
      ) {
        return make_ireturn();
      } else if(CodeUtils.TYPE_LONG.equals(clazz)) {
        return make_lreturn();
      } else if(CodeUtils.TYPE_FLOAT.equals(clazz)) {
        return make_freturn();
      } else if(CodeUtils.TYPE_DOUBLE.equals(clazz)) {
        return make_dreturn();
      }

      return make_areturn();
    }

    /**
     * Make an invoke opcode that is suited to the method.
     *
     * <p>Constructors and private methods are invoked by invokespecial.
     * Static methods are invoked by invokestatic. Methods referred to by
     * interface are invoked by invokeinterface, and methods referred to
     * by classes are invoked by invokeabstract.</p>
     *
     * @param cm  the CodeMethod to invoke
     * @return an invocation Instruction suited to the method
     */
    public static Instruction make_invoke(CodeMethod cm) {
      int modifiers = cm.getModifiers();

      if( (CodeUtils.ACC_STATIC & modifiers) != 0 ) {
        return make_invokestatic(cm);
      } else if( (CodeUtils.ACC_INTERFACE) != 0 ) {
        return make_invokeinterface(cm);
      } else if( (CodeUtils.ACC_PRIVATE) != 0 ) {
        return make_invokespecial(cm);
      } else {
        return make_invokevirtual(cm);
      }
    }
    
    /**
     * Synchronize the processing of an entire block of code on a local
     * variable. The local variable must refer to a Java object. It is safe to
     * replace the value of the local variable within the block.
     *
     * @param lockVar  Label for the local variable to lock on
     * @param code the CodeGenerator that will make the body of the synchronized
     *   block
     * @throws CodeException if locVar holds a primative value
     */
    public static CodeGenerator make_synchronizedBlock(
      LocalVariable lockVar,
      CodeGenerator code
    ) throws CodeException {
      InstructionVector block = new InstructionVector();
      block.add(make_aload(lockVar));
      block.add(make_dup()); // defencive copy incase the local is re-assigned
      block.add(make_monitorenter());
      block.add(code);
      block.add(make_monitorexit());
      
      return block;
    }
    
    /**
     * Synchronize the processing of an entire block of code on the object on
     * the top of the stack.
     *
     * @param lockVar  Label for the local variable to lock on
     * @param code the CodeGenerator that will make the body of the synchronized
     *   block
     */
    public static CodeGenerator make_synchronizedBlock(CodeGenerator code) {
      InstructionVector block = new InstructionVector();
      block.add(make_dup());
      block.add(make_monitorenter());
      block.add(code);
      block.add(make_monitorexit());
      
      return block;
    }
    
    public static CodeGenerator make_markLabel(Label lab) {
      return new MarkLabel(lab);
    }

    /**
     * Make an if Instruction for the opcode and label.
     *
     * <p>Obcodes for IF jump to the label if the condition is true, or execute
     * the next instruction ontherwise. In this case, op will be used to compare
     * the top of the stack and if the condition is true, will cause a jump to
     * the label. Otherwise, the next instruction in the stream will be used.
     * </p>
     *
     * @param op  if opcode
     * @param lab Label to jump to
     * @return an Instruction that will perform an if 
     * @throws IllegalArgumentException if the op code is not an if instruction
     */
    public static Instruction make_if(byte op, Label lab) {
      if(op >= op_ifeq && op <= op_ifle) {
        return new LabelInstruction(op, lab, -1);
      } else if(op >= op_if_icmpeq && op <= op_if_acmpne) {
        return new LabelInstruction(op, lab, -2);
      } else {
        throw new IllegalArgumentException("Opcode must be an if. " + op);
      }
    }
    
    public static ParametricCodeGenerator make_newraray(final ParametricType type) {
      return new ParametricCodeGenerator() {
        public ParametricType getType() { return type; }
        public int stackDepth() { return 0; }
        public int stackDelta() { return 0; }
        public void writeCode(CodeContext ctx) throws CodeException {
          CodeClass cc = ctx.resolveParametricType(type);
          make_newarray(cc).writeCode(ctx);
        }
      };
    }
    
    /**
     * Load an element of a parametric type to an array.
     *
     * <p>This is the parametric version of the make_&lt;x&gt;aload() factory
     * methods.</p>
     *
     * @param type  the ParametricType giving the element type of the array
     * @return a ParametricCodeGenerator that will load the correct type.
     */
    public static ParametricCodeGenerator make_array_load(final ParametricType type) {
      return new ParametricCodeGenerator() {
        public ParametricType getType() { return type; }
        public int stackDepth() { return 0; }
        public int stackDelta() { return -1; }
        public void writeCode(CodeContext ctx) throws CodeException {
          CodeClass cc = ctx.resolveParametricType(type);
          if(!cc.isPrimitive()) {
            ctx.writeByte(op_aaload);
          } else if(cc == CodeUtils.TYPE_INT) {
            ctx.writeByte(op_iaload);
          } else if(cc == CodeUtils.TYPE_LONG) {
            ctx.writeByte(op_daload);
          } else if(cc == CodeUtils.TYPE_FLOAT) {
            ctx.writeByte(op_faload);
          } else if(cc == CodeUtils.TYPE_DOUBLE) {
            ctx.writeByte(op_daload);
          } else if(cc == CodeUtils.TYPE_BYTE) {
            ctx.writeByte(op_baload);
          } else if(cc == CodeUtils.TYPE_CHAR) {
            ctx.writeByte(op_caload);
          } else if(cc == CodeUtils.TYPE_SHORT) {
            ctx.writeByte(op_saload);
          } else {
            throw new CodeException("Confused. Don't recognize type: " + cc);
          }
        }
      };
    }

    /**
     * Store an element of a parametric type to an array.
     *
     * <p>This is the parametric version of the make_&lt;x&gt;astore() factory
     * methods.</p>
     *
     * @param type  the ParametricType giving the element type of the array
     * @return a ParametricCodeGenerator that will store the correct type.
     */
    public static ParametricCodeGenerator make_arrayStore(final ParametricType type) {
      return new ParametricCodeGenerator() {
        public ParametricType getType() { return type; }
        public int stackDepth() { return 0; }
        public int stackDelta() { return -3; }
        public void writeCode(CodeContext ctx) throws CodeException {
          CodeClass cc = ctx.resolveParametricType(type);
          if(!cc.isPrimitive()) {
            ctx.writeByte(op_aastore);
          } else if(cc == CodeUtils.TYPE_INT) {
            ctx.writeByte(op_iastore);
          } else if(cc == CodeUtils.TYPE_LONG) {
            ctx.writeByte(op_dastore);
          } else if(cc == CodeUtils.TYPE_FLOAT) {
            ctx.writeByte(op_fastore);
          } else if(cc == CodeUtils.TYPE_DOUBLE) {
            ctx.writeByte(op_dastore);
          } else if(cc == CodeUtils.TYPE_BYTE) {
            ctx.writeByte(op_bastore);
          } else if(cc == CodeUtils.TYPE_CHAR) {
            ctx.writeByte(op_castore);
          } else if(cc == CodeUtils.TYPE_SHORT) {
            ctx.writeByte(op_sastore);
          } else {
            throw new CodeException("Confused. Don't recognize type: " + cc);
          }
        }
      };
    }
    
    /**
     * Load an item of a parametric type from a local variable.
     *
     * <p>This is the parametric version of the make_&lt;x&gt;load() factory
     * methods.</p>
     *
     * @param type  the ParametricType giving the type to load
     * @param lv  the LocalVariable to load from
     * @return a ParametricCodeGenerator that will load the correct type.
     */
    public ParametricCodeGenerator make_load(
      final ParametricType type,
      final LocalVariable lv
    ) throws CodeException {
      return new ParametricCodeGenerator() {
        public ParametricType getType() { return type; }
        public int stackDepth() { return 1; }
        public int stackDelta() { return 1; }
        public void writeCode(CodeContext ctx) throws CodeException {
          CodeClass cc = ctx.resolveParametricType(type);
          Instruction ins;
          if(!cc.isPrimitive()) {
            ins = make_aload(lv);
          } else if(cc == CodeUtils.TYPE_INT) {
            ins = make_iload(lv);
          } else if(cc == CodeUtils.TYPE_LONG) {
            ins = make_lload(lv);
          } else if(cc == CodeUtils.TYPE_FLOAT) {
            ins = make_fload(lv);
          } else if(cc == CodeUtils.TYPE_DOUBLE) {
            ins = make_dload(lv);
          } else if(cc == CodeUtils.TYPE_BYTE) {
            ins = make_iload(lv);
          } else if(cc == CodeUtils.TYPE_CHAR) {
            ins = make_iload(lv);
          } else if(cc == CodeUtils.TYPE_SHORT) {
            ins = make_iload(lv);
          } else {
            throw new CodeException("Confused. Don't recognize type: " + cc);
          }
          
          ins.writeCode(ctx);
        }
      };
    }
    
    /**
     * Store an item of a parametric type to a local variable.
     *
     * <p>This is the parametric version of the make_&lt;x&gt;save() factory
     * methods.</p>
     *
     * @param type  the ParametricType giving the type to save
     * @param lv  the LocalVariable to load from
     * @return a ParametricCodeGenerator that will save the correct type.
     */
    public ParametricCodeGenerator make_save(
      final ParametricType type,
      final LocalVariable lv
    ) throws CodeException {
      return new ParametricCodeGenerator() {
        public ParametricType getType() { return type; }
        public int stackDepth() { return 0; }
        public int stackDelta() { return -1; }
        public void writeCode(CodeContext ctx) throws CodeException {
          CodeClass cc = ctx.resolveParametricType(type);
          Instruction ins;
          if(!cc.isPrimitive()) {
            ins = make_astore(lv);
          } else if(cc == CodeUtils.TYPE_INT) {
            ins = make_istore(lv);
          } else if(cc == CodeUtils.TYPE_LONG) {
            ins = make_lstore(lv);
          } else if(cc == CodeUtils.TYPE_FLOAT) {
            ins = make_fstore(lv);
          } else if(cc == CodeUtils.TYPE_DOUBLE) {
            ins = make_dstore(lv);
          } else if(cc == CodeUtils.TYPE_BYTE) {
            ins = make_istore(lv);
          } else if(cc == CodeUtils.TYPE_CHAR) {
            ins = make_istore(lv);
          } else if(cc == CodeUtils.TYPE_SHORT) {
            ins = make_istore(lv);
          } else {
            throw new CodeException("Confused. Don't recognize type: " + cc);
          }
          
          ins.writeCode(ctx);
        }
      };
    }
    
    /**
     * Make a return statement for the parametric type.
     *
     * @param type  the ParametricType to return
     */
    public ParametricCodeGenerator make_return(final ParametricType type) {
      return new ParametricCodeGenerator() {
        public ParametricType getType() { return type; }
        public int stackDepth() { return 0; }
        public int stackDelta() { return -1; }
        public void writeCode(CodeContext ctx) throws CodeException {
          CodeClass cc = ctx.resolveParametricType(type);
          if(!cc.isPrimitive()) {
            ctx.writeByte(op_areturn);
          } else if(cc == CodeUtils.TYPE_INT) {
            ctx.writeByte(op_ireturn);
          } else if(cc == CodeUtils.TYPE_LONG) {
            ctx.writeByte(op_lreturn);
          } else if(cc == CodeUtils.TYPE_FLOAT) {
            ctx.writeByte(op_freturn);
          } else if(cc == CodeUtils.TYPE_DOUBLE) {
            ctx.writeByte(op_dreturn);
          } else if(cc == CodeUtils.TYPE_BYTE) {
            ctx.writeByte(op_ireturn);
          } else if(cc == CodeUtils.TYPE_CHAR) {
            ctx.writeByte(op_ireturn);
          } else if(cc == CodeUtils.TYPE_SHORT) {
            ctx.writeByte(op_ireturn);
          } else {
            throw new CodeException("Confused. Don't recognize type: " + cc);
          }
        }
      };
    }
    
    public PParametricCodeGenerator make_cast(
      final ParametricType from,
      final ParametricType to
    ) throws CodeException {
      if(from.isObject()) {
        throw new CodeException("Can not cast from non-primative type: " + from);
      }
      
      if(to.isObject()) {
        throw new CodeException("Can not cast to non-primative type: " + to);
      }
      
      return new PParametricCodeGenerator() {
        public ParametricType getType1() { return from; }
        public ParametricType getType2() { return to; }
        public int stackDepth() { return 0; }
        public int stackDelta() { return 0; }
        public void writeCode(CodeContext ctx) throws CodeException {
          if(from == to) {
            return;
          }
          
          CodeClass fromc = ctx.resolveParametricType(from);
          CodeClass toc = ctx.resolveParametricType(to);
          
          if(!fromc.isPrimitive()) {
            throw new CodeException("Can't cast from non-primitive type: " + fromc);
          }
          
          if(!toc.isPrimitive()) {
            throw new CodeException("Can't cast to non-primitive type: " + toc);
          }
          
          if(fromc == CodeUtils.TYPE_DOUBLE) {
            if(toc == CodeUtils.TYPE_FLOAT) {
              ctx.writeByte(op_d2f);
            } else if(toc == CodeUtils.TYPE_LONG) {
              ctx.writeByte(op_d2l);
            } else {
              ctx.writeByte(op_d2i);
            }
          } else if(fromc == CodeUtils.TYPE_LONG) {
            if(toc == CodeUtils.TYPE_FLOAT) {
              ctx.writeByte(op_l2f);
            } else if(toc == CodeUtils.TYPE_DOUBLE) {
              ctx.writeByte(op_l2d);
            } else {
              ctx.writeByte(op_d2i);
            }
          } else { // something equivalent to integer
            if(toc == CodeUtils.TYPE_FLOAT) {
              ctx.writeByte(op_i2f);
            } else if(toc == CodeUtils.TYPE_DOUBLE) {
              ctx.writeByte(op_i2d);
            } else if(toc == CodeUtils.TYPE_LONG) {
              ctx.writeByte(op_i2l);
            } else if(toc == CodeUtils.TYPE_BYTE) {
              ctx.writeByte(op_i2b);
            } else if(toc == CodeUtils.TYPE_CHAR) {
              ctx.writeByte(op_i2c);
            } else if(toc == CodeUtils.TYPE_SHORT) {
              ctx.writeByte(op_i2s);
            }
          }
        }
      };
    }
    
    public ParametricCodeGenerator make_add(final ParametricType type)
    throws CodeException {
      if(type.isObject()) {
        throw new CodeException("Can't add non-primitive type: " + type);
      }
      
      return new ParametricCodeGenerator() {
        public ParametricType getType() { return type; }
        public int stackDelta() { return -1; }
        public int stackDepth() { return 0; }
        public void writeCode(CodeContext ctx) throws CodeException {
          CodeClass cc = ctx.resolveParametricType(type);
          
          if(!cc.isPrimitive()) {
            throw new CodeException(
              "We can only add primitive types: " +
              type + " : " + cc );
          }
          
          if(cc == CodeUtils.TYPE_DOUBLE) {
            ctx.writeByte(op_dadd);
          } else if(cc == CodeUtils.TYPE_LONG) {
            ctx.writeByte(op_ladd);
          } else if(cc == CodeUtils.TYPE_FLOAT) {
            ctx.writeByte(op_fadd);
          } else {
            ctx.writeByte(op_iadd);
          }
        }
      };
    }
    
    public ParametricCodeGenerator make_sub(final ParametricType type)
    throws CodeException {
      if(type.isObject()) {
        throw new CodeException("Can't sub non-primitive type: " + type);
      }
      
      return new ParametricCodeGenerator() {
        public ParametricType getType() { return type; }
        public int stackDelta() { return -1; }
        public int stackDepth() { return 0; }
        public void writeCode(CodeContext ctx) throws CodeException {
          CodeClass cc = ctx.resolveParametricType(type);
          
          if(!cc.isPrimitive()) {
            throw new CodeException(
              "We can only sub primitive types: " +
              type + " : " + cc );
          }
          
          if(cc == CodeUtils.TYPE_DOUBLE) {
            ctx.writeByte(op_dsub);
          } else if(cc == CodeUtils.TYPE_LONG) {
            ctx.writeByte(op_lsub);
          } else if(cc == CodeUtils.TYPE_FLOAT) {
            ctx.writeByte(op_fsub);
          } else {
            ctx.writeByte(op_isub);
          }
        }
      };
    }

    
    public ParametricCodeGenerator make_mul(final ParametricType type)
    throws CodeException {
      if(type.isObject()) {
        throw new CodeException("Can't mul non-primitive type: " + type);
      }
      
      return new ParametricCodeGenerator() {
        public ParametricType getType() { return type; }
        public int stackDelta() { return -1; }
        public int stackDepth() { return 0; }
        public void writeCode(CodeContext ctx) throws CodeException {
          CodeClass cc = ctx.resolveParametricType(type);
          
          if(!cc.isPrimitive()) {
            throw new CodeException(
              "We can only mul primitive types: " +
              type + " : " + cc );
          }
          
          if(cc == CodeUtils.TYPE_DOUBLE) {
            ctx.writeByte(op_dmul);
          } else if(cc == CodeUtils.TYPE_LONG) {
            ctx.writeByte(op_lmul);
          } else if(cc == CodeUtils.TYPE_FLOAT) {
            ctx.writeByte(op_fmul);
          } else {
            ctx.writeByte(op_imul);
          }
        }
      };
    }
    
    
    public ParametricCodeGenerator make_div(final ParametricType type)
    throws CodeException {
      if(type.isObject()) {
        throw new CodeException("Can't div non-primitive type: " + type);
      }
      
      return new ParametricCodeGenerator() {
        public ParametricType getType() { return type; }
        public int stackDelta() { return -1; }
        public int stackDepth() { return 0; }
        public void writeCode(CodeContext ctx) throws CodeException {
          CodeClass cc = ctx.resolveParametricType(type);
          
          if(!cc.isPrimitive()) {
            throw new CodeException(
              "We can only div primitive types: " +
              type + " : " + cc );
          }
          
          if(cc == CodeUtils.TYPE_DOUBLE) {
            ctx.writeByte(op_ddiv);
          } else if(cc == CodeUtils.TYPE_LONG) {
            ctx.writeByte(op_ldiv);
          } else if(cc == CodeUtils.TYPE_FLOAT) {
            ctx.writeByte(op_fdiv);
          } else {
            ctx.writeByte(op_idiv);
          }
        }
      };
    }
    
    public ParametricCodeGenerator make_rem(final ParametricType type)
    throws CodeException {
      if(type.isObject()) {
        throw new CodeException("Can't rem non-primitive type: " + type);
      }
      
      return new ParametricCodeGenerator() {
        public ParametricType getType() { return type; }
        public int stackDelta() { return -1; }
        public int stackDepth() { return 0; }
        public void writeCode(CodeContext ctx) throws CodeException {
          CodeClass cc = ctx.resolveParametricType(type);
          
          if(!cc.isPrimitive()) {
            throw new CodeException(
              "We can only rem primitive types: " +
              type + " : " + cc );
          }
          
          if(cc == CodeUtils.TYPE_DOUBLE) {
            ctx.writeByte(op_drem);
          } else if(cc == CodeUtils.TYPE_LONG) {
            ctx.writeByte(op_lrem);
          } else if(cc == CodeUtils.TYPE_FLOAT) {
            ctx.writeByte(op_frem);
          } else {
            ctx.writeByte(op_irem);
          }
        }
      };
    }
    
    public ParametricCodeGenerator make_neg(final ParametricType type)
    throws CodeException {
      if(type.isObject()) {
        throw new CodeException("Can't add non-primitive type: " + type);
      }
      
      return new ParametricCodeGenerator() {
        public ParametricType getType() { return type; }
        public int stackDelta() { return 0; }
        public int stackDepth() { return 0; }
        public void writeCode(CodeContext ctx) throws CodeException {
          CodeClass cc = ctx.resolveParametricType(type);
          
          if(!cc.isPrimitive()) {
            throw new CodeException(
              "We can only add primitive types: " +
              type + " : " + cc );
          }
          
          if(cc == CodeUtils.TYPE_DOUBLE) {
            ctx.writeByte(op_dneg);
          } else if(cc == CodeUtils.TYPE_LONG) {
            ctx.writeByte(op_lneg);
          } else if(cc == CodeUtils.TYPE_FLOAT) {
            ctx.writeByte(op_fneg);
          } else {
            ctx.writeByte(op_ineg);
          }
        }
      };
    }

    public ParametricCodeGenerator make_shiftLeft(final ParametricType type)
    throws CodeException {
      if(type.isObject()) {
        throw new CodeException("Can't shift non-primitive type: " + type);
      }
      return new ParametricCodeGenerator() {
        public ParametricType getType() { return type; }
        public int stackDelta() { return -1; }
        public int stackDepth() { return 0; }
        public void writeCode(CodeContext ctx) throws CodeException {
          CodeClass cc = ctx.resolveParametricType(type);
          
          if(!cc.isPrimitive()) {
            throw new CodeException("Can't shift non-primitive type: " + cc);
          }
          
          if(CodeUtils.isFloatType(cc)) {
            throw new CodeException("Can't shift floating point type: " + cc);
          }
          
          if(cc == CodeUtils.TYPE_LONG) {
            ctx.writeByte(op_lshl);
          } else  {
            ctx.writeByte(op_ishl);
          }
        }
      };
    }
    
    public ParametricCodeGenerator make_shiftRight(final ParametricType type)
    throws CodeException {
      if(type.isObject()) {
        throw new CodeException("Can't shift non-primitive type: " + type);
      }
      return new ParametricCodeGenerator() {
        public ParametricType getType() { return type; }
        public int stackDelta() { return -1; }
        public int stackDepth() { return 0; }
        public void writeCode(CodeContext ctx) throws CodeException {
          CodeClass cc = ctx.resolveParametricType(type);
          
          if(!cc.isPrimitive()) {
            throw new CodeException("Can't shift non-primitive type: " + cc);
          }
          
          if(CodeUtils.isFloatType(cc)) {
            throw new CodeException("Can't shift floating point type: " + cc);
          }
          
          if(cc == CodeUtils.TYPE_LONG) {
            ctx.writeByte(op_lshr);
          } else {
            ctx.writeByte(op_ishr);
          }
        }
      };
    }
    
    public ParametricCodeGenerator make_shiftRightLogical(final ParametricType type)
    throws CodeException {
      if(type.isObject()) {
        throw new CodeException("Can't shift non-primitive type: " + type);
      }
      return new ParametricCodeGenerator() {
        public ParametricType getType() { return type; }
        public int stackDelta() { return -1; }
        public int stackDepth() { return 0; }
        public void writeCode(CodeContext ctx) throws CodeException {
          CodeClass cc = ctx.resolveParametricType(type);
          
          if(!cc.isPrimitive()) {
            throw new CodeException("Can't shift non-primitive type: " + cc);
          }
          
          if(CodeUtils.isFloatType(cc)) {
            throw new CodeException("Can't shift floating point type: " + cc);
          }
          
          if(cc == CodeUtils.TYPE_LONG) {
            ctx.writeByte(op_lushr);
          } else {
            ctx.writeByte(op_iushr);
          }
        }
      };
    }
    
    public ParametricCodeGenerator make_and(final ParametricType type)
    throws CodeException {
      if(type.isObject()) {
        throw new CodeException("Can't and non-primitive type: " + type);
      }
      return new ParametricCodeGenerator() {
        public ParametricType getType() { return type; }
        public int stackDelta() { return -1; }
        public int stackDepth() { return 0; }
        public void writeCode(CodeContext ctx) throws CodeException {
          CodeClass cc = ctx.resolveParametricType(type);
          
          if(!cc.isPrimitive()) {
            throw new CodeException("Can't and non-primitive type: " + cc);
          }
          
          if(CodeUtils.isFloatType(cc)) {
            throw new CodeException("Can't and floating point type: " + cc);
          }
          
          
          if(cc == CodeUtils.TYPE_LONG) {
            ctx.writeByte(op_land);
          } else {
            ctx.writeByte(op_iand);
          }
        }
      };
    }
    
    public ParametricCodeGenerator make_or(final ParametricType type)
    throws CodeException {
      if(type.isObject()) {
        throw new CodeException("Can't and non-primitive type: " + type);
      }
      return new ParametricCodeGenerator() {
        public ParametricType getType() { return type; }
        public int stackDelta() { return -1; }
        public int stackDepth() { return 0; }
        public void writeCode(CodeContext ctx) throws CodeException {
          CodeClass cc = ctx.resolveParametricType(type);
          
          if(!cc.isPrimitive()) {
            throw new CodeException("Can't or non-primitive type: " + cc);
          }
          
          if(CodeUtils.isFloatType(cc)) {
            throw new CodeException("Can't or floating point type: " + cc);
          }
          
          
          if(cc == CodeUtils.TYPE_LONG) {
            ctx.writeByte(op_lor);
          } else {
            ctx.writeByte(op_ior);
          }
        }
      };
    }
    
    public ParametricCodeGenerator make_xor(final ParametricType type)
    throws CodeException {
      if(type.isObject()) {
        throw new CodeException("Can't and non-primitive type: " + type);
      }
      return new ParametricCodeGenerator() {
        public ParametricType getType() { return type; }
        public int stackDelta() { return -1; }
        public int stackDepth() { return 0; }
        public void writeCode(CodeContext ctx) throws CodeException {
          CodeClass cc = ctx.resolveParametricType(type);
          
          if(!cc.isPrimitive()) {
            throw new CodeException("Can't xor non-primitive type: " + cc);
          }

          if(CodeUtils.isFloatType(cc)) {
            throw new CodeException("Can't xor floating point type: " + cc);
          }
          
          
          if(cc == CodeUtils.TYPE_LONG) {
            ctx.writeByte(op_lxor);
          } else  {
            ctx.writeByte(op_ixor);
          }
        }
      };
    }
    
    public ParametricCodeGenerator make_cmpg(final ParametricType type)
    throws CodeException {
      if(type.isObject()) {
        throw new CodeException("Can't cmpg non-primitive type: " + type);
      }
      
      return new ParametricCodeGenerator() {
        public ParametricType getType() { return type; }
        public int stackDelta() { return -1; }
        public int stackDepth() { return 0; }
        public void writeCode(CodeContext ctx) throws CodeException {
          CodeClass cc = ctx.resolveParametricType(type);
          
          if(!CodeUtils.isFloatType(cc)) {
            throw new CodeException("Can only cmpg floating point types: " + cc);
          }
          
          if(cc == CodeUtils.TYPE_DOUBLE) {
            ctx.writeByte(op_dcmpg);
          } else  {
            ctx.writeByte(op_fcmpg);
          }
        }
      };
    }
    
    public ParametricCodeGenerator make_cmpl(final ParametricType type)
    throws CodeException {
      if(type.isObject()) {
        throw new CodeException("Can't cmpl non-primitive type: " + type);
      }
      
      return new ParametricCodeGenerator() {
        public ParametricType getType() { return type; }
        public int stackDelta() { return -1; }
        public int stackDepth() { return 0; }
        public void writeCode(CodeContext ctx) throws CodeException {
          CodeClass cc = ctx.resolveParametricType(type);
          
          if(!CodeUtils.isFloatType(cc)) {
            throw new CodeException("Can only cmpl floating point types: " + cc);
          }
          
          if(cc == CodeUtils.TYPE_DOUBLE) {
            ctx.writeByte(op_dcmpl);
          } else  {
            ctx.writeByte(op_fcmpl);
          }
        }
      };
    }
    
    public ParametricCodeGenerator make_dup(final ParametricType type)
    throws CodeException {
      return new ParametricCodeGenerator() {
        public ParametricType getType() { return type; }
        public int stackDelta() { return 1; }
        public int stackDepth() { return 1; }
        public void writeCode(CodeContext ctx) throws CodeException {
          CodeClass cc = ctx.resolveParametricType(type);
          
          if(CodeUtils.wordsForType(cc) == 2) {
            ctx.writeByte(op_dup2);
          } else  {
            ctx.writeByte(op_dup);
          }
        }
      };
    }
    
    //
    // Java opcodes
    //

    public final static byte op_nop = 0;
    public final static byte op_aconst_null = 1;
    public final static byte op_iconst_m1 = 2;
    public final static byte op_iconst_0 = 3;
    public final static byte op_iconst_1 = 4;
    public final static byte op_iconst_2 = 5;
    public final static byte op_iconst_3 = 6;
    public final static byte op_iconst_4 = 7;
    public final static byte op_iconst_5 = 8;
    public final static byte op_lconst_0 = 9;
    public final static byte op_lconst_1 = 10;
    public final static byte op_fconst_0 = 11;
    public final static byte op_fconst_1 = 12;
    public final static byte op_fconst_2 = 13;
    public final static byte op_dconst_0 = 14;
    public final static byte op_dconst_1 = 15;
    public final static byte op_bipush = 16;
    public final static byte op_sipush = 17;
    public final static byte op_ldc = 18;
    public final static byte op_ldc_w = 19;
    public final static byte op_ldc2_w = 20;
    public final static byte op_iload = 21;
    public final static byte op_lload = 22;
    public final static byte op_fload = 23;
    public final static byte op_dload = 24;
    public final static byte op_aload = 25;
    public final static byte op_iload_0 = 26;
    public final static byte op_iload_1 = 27;
    public final static byte op_iload_2 = 28;
    public final static byte op_iload_3 = 29;
    public final static byte op_lload_0 = 30;
    public final static byte op_lload_1 = 31;
    public final static byte op_lload_2 = 32;
    public final static byte op_lload_3 = 33;
    public final static byte op_fload_0 = 34;
    public final static byte op_fload_1 = 35;
    public final static byte op_fload_2 = 36;
    public final static byte op_fload_3 = 37;
    public final static byte op_dload_0 = 38;
    public final static byte op_dload_1 = 39;
    public final static byte op_dload_2 = 40;
    public final static byte op_dload_3 = 41;
    public final static byte op_aload_0 = 42;
    public final static byte op_aload_1 = 43;
    public final static byte op_aload_2 = 44;
    public final static byte op_aload_3 = 45;
    public final static byte op_iaload = 46;
    public final static byte op_laload = 47;
    public final static byte op_faload = 48;
    public final static byte op_daload = 49;
    public final static byte op_aaload = 50;
    public final static byte op_baload = 51;
    public final static byte op_caload = 52;
    public final static byte op_saload = 53;
    public final static byte op_istore = 54;
    public final static byte op_lstore = 55;
    public final static byte op_fstore = 56;
    public final static byte op_dstore = 57;
    public final static byte op_astore = 58;
    public final static byte op_istore_0 = 59;
    public final static byte op_istore_1 = 60;
    public final static byte op_istore_2 = 61;
    public final static byte op_istore_3 = 62;
    public final static byte op_lstore_0 = 63;
    public final static byte op_lstore_1 = 64;
    public final static byte op_lstore_2 = 65;
    public final static byte op_lstore_3 = 66;
    public final static byte op_fstore_0 = 67;
    public final static byte op_fstore_1 = 68;
    public final static byte op_fstore_2 = 69;
    public final static byte op_fstore_3 = 70;
    public final static byte op_dstore_0 = 71;
    public final static byte op_dstore_1 = 72;
    public final static byte op_dstore_2 = 73;
    public final static byte op_dstore_3 = 74;
    public final static byte op_astore_0 = 75;
    public final static byte op_astore_1 = 76;
    public final static byte op_astore_2 = 77;
    public final static byte op_astore_3 = 78;
    public final static byte op_iastore = 79;
    public final static byte op_lastore = 80;
    public final static byte op_fastore = 81;
    public final static byte op_dastore = 82;
    public final static byte op_aastore = 83;
    public final static byte op_bastore = 84;
    public final static byte op_castore = 85;
    public final static byte op_sastore = 86;
    public final static byte op_pop = 87;
    public final static byte op_pop2 = 88;
    public final static byte op_dup = 89;
    public final static byte op_dup_x1 = 90;
    public final static byte op_dup_x2 = 91;
    public final static byte op_dup2 = 92;
    public final static byte op_dup2_x1 = 93;
    public final static byte op_dup2_x2 = 94;
    public final static byte op_swap = 95;
    public final static byte op_iadd = 96;
    public final static byte op_ladd = 97;
    public final static byte op_fadd = 98;
    public final static byte op_dadd = 99;
    public final static byte op_isub = 100;
    public final static byte op_lsub = 101;
    public final static byte op_fsub = 102;
    public final static byte op_dsub = 103;
    public final static byte op_imul = 104;
    public final static byte op_lmul = 105;
    public final static byte op_fmul = 106;
    public final static byte op_dmul = 107;
    public final static byte op_idiv = 108;
    public final static byte op_ldiv = 109;
    public final static byte op_fdiv = 110;
    public final static byte op_ddiv = 111;
    public final static byte op_irem = 112;
    public final static byte op_lrem = 113;
    public final static byte op_frem = 114;
    public final static byte op_drem = 115;
    public final static byte op_ineg = 116;
    public final static byte op_lneg = 117;
    public final static byte op_fneg = 118;
    public final static byte op_dneg = 119;
    public final static byte op_ishl = 120;
    public final static byte op_lshl = 121;
    public final static byte op_ishr = 122;
    public final static byte op_lshr = 123;
    public final static byte op_iushr = 124;
    public final static byte op_lushr = 125;
    public final static byte op_iand = 126;
    public final static byte op_land = 127;
    public final static byte op_ior = (byte) 128;
    public final static byte op_lor = (byte) 129;
    public final static byte op_ixor = (byte) 130;
    public final static byte op_lxor = (byte) 131;
    public final static byte op_iinc = (byte) 132;
    public final static byte op_i2l = (byte) 133;
    public final static byte op_i2f = (byte) 134;
    public final static byte op_i2d = (byte) 135;
    public final static byte op_l2i = (byte) 136;
    public final static byte op_l2f = (byte) 137;
    public final static byte op_l2d = (byte) 138;
    public final static byte op_f2i = (byte) 139;
    public final static byte op_f2l = (byte) 140;
    public final static byte op_f2d = (byte) 141;
    public final static byte op_d2i = (byte) 142;
    public final static byte op_d2l = (byte) 143;
    public final static byte op_d2f = (byte) 144;
    public final static byte op_i2b = (byte) 145;
    public final static byte op_i2c = (byte) 146;
    public final static byte op_i2s = (byte) 147;
    public final static byte op_lcmp = (byte) 148;
    public final static byte op_fcmpl = (byte) 149;
    public final static byte op_fcmpg = (byte) 150;
    public final static byte op_dcmpl = (byte) 151;
    public final static byte op_dcmpg = (byte) 152;
    public final static byte op_ifeq = (byte) 153;
    public final static byte op_ifne = (byte) 154;
    public final static byte op_iflt = (byte) 155;
    public final static byte op_ifge = (byte) 156;
    public final static byte op_ifgt = (byte) 157;
    public final static byte op_ifle = (byte) 158;
    public final static byte op_if_icmpeq = (byte) 159;
    public final static byte op_if_icmpne = (byte) 160;
    public final static byte op_if_icmplt = (byte) 161;
    public final static byte op_if_icmpge = (byte) 162;
    public final static byte op_if_icmpgt = (byte) 163;
    public final static byte op_if_icmple = (byte) 164;
    public final static byte op_if_acmpeq = (byte) 165;
    public final static byte op_if_acmpne = (byte) 166;
    public final static byte op_goto = (byte) 167;
    public final static byte op_jsr = (byte) 168;
    public final static byte op_ret = (byte) 169;
    public final static byte op_tableswitch = (byte) 170;
    public final static byte op_lookupswitch = (byte) 171;
    public final static byte op_ireturn = (byte) 172;
    public final static byte op_lreturn = (byte) 173;
    public final static byte op_freturn = (byte) 174;
    public final static byte op_dreturn = (byte) 175;
    public final static byte op_areturn = (byte) 176;
    public final static byte op_return = (byte) 177;
    public final static byte op_getstatic = (byte) 178;
    public final static byte op_putstatic = (byte) 179;
    public final static byte op_getfield = (byte) 180;
    public final static byte op_putfield = (byte) 181;
    public final static byte op_invokevirtual = (byte) 182;
    public final static byte op_invokespecial = (byte) 183;
    public final static byte op_invokestatic = (byte) 184;
    public final static byte op_invokeinterface = (byte) 185;
    public final static byte op_new = (byte) 187;
    public final static byte op_newarray = (byte) 188;
    public final static byte op_anewarray = (byte) 189;
    public final static byte op_arraylength = (byte) 190;
    public final static byte op_athrow = (byte) 191;
    public final static byte op_checkcast = (byte) 192;
    public final static byte op_instanceof = (byte) 193;
    public final static byte op_monitorenter = (byte) 194;
    public final static byte op_monitorexit = (byte) 195;
    public final static byte op_wide = (byte) 196;
    public final static byte op_multianewarray = (byte) 197;
    public final static byte op_ifnull = (byte) 198;
    public final static byte op_ifnonnull = (byte) 199;
    public final static byte op_goto_w = (byte) 200;
    public final static byte op_jsr_w = (byte) 201;
    public final static byte op_breakpoint = (byte) 202;
    public final static byte op_ldc_quick = (byte) 203;
    public final static byte op_ldc_w_quick = (byte) 204;
    public final static byte op_ldc2_w_quick = (byte) 205;
    public final static byte op_getfield_quick = (byte) 206;
    public final static byte op_putfield_quick = (byte) 207;
    public final static byte op_getfield2_quick = (byte) 208;
    public final static byte op_putfield2_quick = (byte) 209;
    public final static byte op_getstatic_quick = (byte) 210;
    public final static byte op_putstatic_quick = (byte) 211;
    public final static byte op_getstatic2_quick = (byte) 212;
    public final static byte op_putstatic2_quick = (byte) 213;
    public final static byte op_invokevirtual_quick = (byte) 214;
    public final static byte op_invokenonvirtual_quick = (byte) 215;
    public final static byte op_invokesuper_quick = (byte) 216;
    public final static byte op_invokestatic_quick = (byte) 217;
    public final static byte op_invokeinterface_quick = (byte) 218;
    public final static byte op_invokevirtualobject_quick = (byte) 219;
    public final static byte op_new_quick = (byte) 221;
    public final static byte op_anewarray_quick = (byte) 222;
    public final static byte op_multianewarray_quick = (byte) 223;
    public final static byte op_checkcast_quick = (byte) 224;
    public final static byte op_instanceof_quick = (byte) 225;
    public final static byte op_invokevirtual_quick_w = (byte) 226;
    public final static byte op_getfield_quick_w = (byte) 227;
    public final static byte op_putfield_quick_w = (byte) 228;
    public final static byte op_impdep1 = (byte) 254;
    public final static byte op_impdep2 = (byte) 255;
}
