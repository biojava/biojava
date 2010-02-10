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


package org.biojava.bio.dp.twohead;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.BackPointer;
import org.biojava.bio.dp.DP;
import org.biojava.bio.dp.EmissionState;
import org.biojava.bio.dp.IllegalTransitionException;
import org.biojava.bio.dp.MagicalState;
import org.biojava.bio.dp.MarkovModel;
import org.biojava.bio.dp.ScoreType;
import org.biojava.bio.dp.State;
import org.biojava.bio.symbol.AlphabetIndex;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.ClassTools;
import org.biojava.utils.bytecode.ByteCode;
import org.biojava.utils.bytecode.CodeClass;
import org.biojava.utils.bytecode.CodeException;
import org.biojava.utils.bytecode.CodeField;
import org.biojava.utils.bytecode.CodeGenerator;
import org.biojava.utils.bytecode.CodeMethod;
import org.biojava.utils.bytecode.CodeUtils;
import org.biojava.utils.bytecode.GeneratedClassLoader;
import org.biojava.utils.bytecode.GeneratedCodeClass;
import org.biojava.utils.bytecode.GeneratedCodeMethod;
import org.biojava.utils.bytecode.IfExpression;
import org.biojava.utils.bytecode.InstructionVector;
import org.biojava.utils.bytecode.IntrospectedCodeClass;
import org.biojava.utils.bytecode.LocalVariable;

/**
 * This is an implementation of CellCalculatorFactoryMaker that compiles the
 * HMM object down to Java byte-code that is equivalent in behaviour to the
 * interpreter.
 *
 * @author Matthew Pocock
 * @author Greg Cox
 * @since 1.1
 */

public class DPCompiler implements CellCalculatorFactoryMaker {
  //private final boolean debug = false;
  private final boolean debug = true;

  private final boolean dumpToDisk;

  private final GeneratedClassLoader classLoader;

  public DPCompiler(boolean dumpToDisk) {
    this.dumpToDisk = dumpToDisk;
    classLoader = new GeneratedClassLoader(ClassTools.getClassLoader(this));
  }

  public CellCalculatorFactory make(DP dp) {
    Class forwardC = generateForardClass(dp);
    Class backwardC = generateBackwardClass(dp);
    Class viterbiC = generateViterbiClass(dp);

    try {
      Constructor forward = forwardC.getConstructor(new Class[] {
        DP.class,
        ScoreType.class
      });
      Constructor backward = backwardC.getConstructor(new Class[] {
        DP.class,
        ScoreType.class
      });
      Constructor viterbi = viterbiC.getConstructor(new Class[] {
        DP.class,
        ScoreType.class,
        BackPointer.class
      });
      return new Factory(dp, forward, backward, viterbi);
    } catch (NoSuchMethodException nsme) {
      throw new BioError("Couldn't find constructor on generated class", nsme);
    }
  }

  public static String makeName(String prefix, MarkovModel model) {
    StringBuffer nameBuffer = new StringBuffer(prefix);

    for(Iterator i = model.stateAlphabet().iterator(); i.hasNext(); ) {
      nameBuffer.append("_");
      try {
        nameBuffer.append(model.transitionsFrom((State) i.next()).size());
      } catch (IllegalSymbolException ise) {
        throw new BioError(
          "Assertion Failure: State dissapeared from model", ise
        );
      }
    }

    return nameBuffer.substring(0);
  }

  public Class generateForardClass(DP dp) {
    // forward recursion is:
    //
    // score[state_i] =
    //   emission * sum_j(prev_score[state_j] * transition[state_j, state_i])
    //
    // where prev_score is found as score[advance_i[0][advance_i[1]]
    // however, in log space (as we use), this becomes:
    //
    // score[state_i] =
    //   emission + log( sum_j (
    //     exp(prev_score[state_j] + transition[state_i, state_j])
    //   ))
    //
    // In practice, for sequences of any length, the sum terms are too
    // near to zero for numerical instablilty not to play a huge part. So, we
    // take out a factor of max_j( prev_score[state_j] ) from the sum and add
    // it to the total.
    //
    // score[state_i] =
    //   emission + log( sum_j (
    //     exp(transition[state_i, state_j] + (prev_score[state_j]- factor) )
    //   )) + factor


    try {
      MarkovModel model = dp.getModel();

      String name = makeName("org.biojava.bio.dp.twohead.Forward", model);
      if(!classLoader.hasGeneratedClass(name)) {
        CodeClass _Object = IntrospectedCodeClass.forClass(Object.class);
        CodeClass _DP = IntrospectedCodeClass.forClass(DP.class);
        CodeClass _ScoreType = IntrospectedCodeClass.forClass(ScoreType.class);
        CodeClass _State_A = IntrospectedCodeClass.forClass(State [].class);
        CodeClass _Cell_A_A = IntrospectedCodeClass.forClass(Cell [][].class);
        CodeClass _CellCalculator = IntrospectedCodeClass.forClass(CellCalculator.class);

        GeneratedCodeClass clazz = new GeneratedCodeClass(
          name,
          _Object,
          new CodeClass[] {_CellCalculator},
          CodeUtils.ACC_PUBLIC | CodeUtils.ACC_STRICT
        );

        State[] states = dp.getStates();
        AlphabetIndex stateIndexer = AlphabetManager.getAlphabetIndex(states);

        CodeField stateF = clazz.createField(
          "states",
          _State_A,
          CodeUtils.ACC_PROTECTED
        );

        // The fields that contain transition scores as double [].
        // This is built so that if two states shair the same transition
        // distribution, they refer to the same transition field. This gives
        // good optimizers a fighting chance to do some hard-core optimizations.
        CodeField[] transitionFields = new CodeField[states.length];
        AlphabetIndex[] transitionIndexers = new AlphabetIndex[states.length];

        // The constructor must load in the transition probabilities.
        // It uses the indexToFieldIndex to ensure that parameters are loaded
        // just once.
        GeneratedCodeMethod __init = clazz.createMethod(
          "<init>",
          CodeUtils.TYPE_VOID,
          new CodeClass[] {_DP, _ScoreType },
          new String[] { "dp", "scoreType" },
          CodeUtils.ACC_PUBLIC
        );

        clazz.setCodeGenerator(__init, createInit(
          model, states,
          clazz, __init,
          transitionIndexers, transitionFields,
          stateF, null
        ));

        GeneratedCodeMethod initialize = clazz.createMethod(
          "initialize",
          CodeUtils.TYPE_VOID,
          new CodeClass[] {_Cell_A_A},
          new String[] { "cells" },
          CodeUtils.ACC_PUBLIC
        );
        clazz.setCodeGenerator(
          initialize,
          createFRecursion(
            true,
            model, states, stateIndexer,
            stateF, transitionFields, transitionIndexers,
            initialize
          )
        );

        GeneratedCodeMethod calcCell = clazz.createMethod(
          "calcCell",
          CodeUtils.TYPE_VOID,
          new CodeClass[] {_Cell_A_A },
          new String[] { "cells" },
          CodeUtils.ACC_PUBLIC
        );
        clazz.setCodeGenerator(
        calcCell,
        createFRecursion(
          false,
          model, states, stateIndexer,
          stateF, transitionFields, transitionIndexers,
          calcCell
        )
      );

      if(dumpToDisk == true) {
        try {
          StringBuffer fName = new StringBuffer(name);
          for(int i = 0; i < fName.length(); i++) {
            if(fName.charAt(i) == '.') {
              fName.setCharAt(i, '/');
            }
          }
          fName.append(".class");

          File dumpFile = new File(fName.substring(0));
          System.out.println("Dumping file: " + fName);
          dumpFile.getParentFile().mkdirs();

          OutputStream out = new FileOutputStream(dumpFile);
          clazz.createCode(out);
          out.close();
        } catch (FileNotFoundException fnfe) {
          throw new BioError("Couldn't dump dp class", fnfe);
        } catch (IOException ioe) {
          throw new BioError("Couldn't dump dp class", ioe);
        }
      }

      classLoader.defineClass(clazz);
    }
     try {
          return classLoader.loadClass(name);
        } catch (Exception e) {
          throw new BioError("Can't find previously generated class for " + name, e);
        }
    } catch (CodeException ce) {
      throw new BioError("Couldn't generate class", ce);
    } catch (NoSuchMethodException nsme) {
      throw new BioError( "Couldn't find method", nsme);
    } catch (NoSuchFieldException nsfe) {
      throw new BioError("Couldn't find field", nsfe);
    } catch (IllegalSymbolException ise) {
      throw new BioError("Couldn't find symbol", ise);
    } catch (BioException be) {
      throw new BioError("Couldn't create indexer", be);
    }
  }

  public Class generateBackwardClass(DP dp) {
    // backward recursion is:
    // score[state_i] =
    //   sum_j(prev_score[state_j] * transition[state_i, state_j] * emission[j])
    // where prev_score is found as score[advance_j[0][advance_j[1]]
    // however, in log space (as we use), this becomes:
    // score[state_i] =
    //   emission + log( sum_j (
    //     exp(prev_score[state_j] + transition[state_j, state_i)
    //   ))
    //
    // In practice, for sequences of any length, the sum terms are too
    // near to zero for numerical instablilty not to play a huge part. So, we
    // take out a factor of max_j( prev_score[state_j] ) from the sum and add
    // it to the total.

    try {
      MarkovModel model = dp.getModel();

      String name = makeName("org.biojava.bio.dp.twohead.Backward", model);
      if(!classLoader.hasGeneratedClass(name)) {
        CodeClass _Object = IntrospectedCodeClass.forClass(Object.class);
        CodeClass _DP = IntrospectedCodeClass.forClass(DP.class);
        CodeClass _ScoreType = IntrospectedCodeClass.forClass(ScoreType.class);
        CodeClass _State_A = IntrospectedCodeClass.forClass(State [].class);
        CodeClass _Cell_A_A = IntrospectedCodeClass.forClass(Cell [][].class);
        CodeClass _CellCalculator = IntrospectedCodeClass.forClass(CellCalculator.class);

        GeneratedCodeClass clazz = new GeneratedCodeClass(
          name,
          _Object,
          new CodeClass[] {_CellCalculator},
          CodeUtils.ACC_PUBLIC | CodeUtils.ACC_STRICT
        );

        State[] states = dp.getStates();
        AlphabetIndex stateIndexer = AlphabetManager.getAlphabetIndex(states);

        CodeField stateF = clazz.createField(
          "states",
          _State_A,
          CodeUtils.ACC_PROTECTED
        );

        // The fields that contain transition scores as double [].
        // This is built so that if two states shair the same transition
        // distribution, they refer to the same transition field. This gives
        // good optimizers a fighting chance to do some hard-core optimizations.
        CodeField[] transitionFields = new CodeField[states.length];
        AlphabetIndex[] transitionIndexers = new AlphabetIndex[states.length];

        // The constructor must load in the transition probabilities.
        // It uses the indexToFieldIndex to ensure that parameters are loaded
        // just once.
        GeneratedCodeMethod __init = clazz.createMethod(
          "<init>",
          CodeUtils.TYPE_VOID,
          new CodeClass[] {_DP, _ScoreType },
          new String[] { "dp", "scoreType" },
          CodeUtils.ACC_PUBLIC
        );

        clazz.setCodeGenerator(__init, createInit(
          model, states,
          clazz, __init,
          transitionIndexers, transitionFields,
          stateF, null
        ));

        GeneratedCodeMethod initialize = clazz.createMethod(
          "initialize",
          CodeUtils.TYPE_VOID,
          new CodeClass[] {_Cell_A_A},
          new String[] { "cells" },
          CodeUtils.ACC_PUBLIC
        );
        clazz.setCodeGenerator(
          initialize,
          createBRecursion(
            true,
            model, states, stateIndexer,
            stateF, transitionFields, transitionIndexers,
            initialize
          )
        );

        GeneratedCodeMethod calcCell = clazz.createMethod(
          "calcCell",
          CodeUtils.TYPE_VOID,
          new CodeClass[] {_Cell_A_A },
          new String[] { "cells" },
          CodeUtils.ACC_PUBLIC
        );
        clazz.setCodeGenerator(
        calcCell,
        createBRecursion(
          false,
          model, states, stateIndexer,
          stateF, transitionFields, transitionIndexers,
          calcCell
        )
      );

      if(dumpToDisk == true) {
        try {
          StringBuffer fName = new StringBuffer(name);
          for(int i = 0; i < fName.length(); i++) {
            if(fName.charAt(i) == '.') {
              fName.setCharAt(i, '/');
            }
          }
          fName.append(".class");

          File dumpFile = new File(fName.substring(0));
          System.out.println("Dumping file: " + fName);
          dumpFile.getParentFile().mkdirs();

          OutputStream out = new FileOutputStream(dumpFile);
          clazz.createCode(out);
          out.close();
        } catch (FileNotFoundException fnfe) {
          throw new BioError( "Couldn't dump dp class", fnfe);
        } catch (IOException ioe) {
          throw new BioError( "Couldn't dump dp class", ioe);
        }
      }

      classLoader.defineClass(clazz);
    }
     try {
          return classLoader.loadClass(name);
        } catch (Exception e) {
          throw new BioError("Can't find previously generated class for " + name, e);
        }
    } catch (CodeException ce) {
      throw new BioError("Couldn't generate class", ce);
    } catch (NoSuchMethodException nsme) {
      throw new BioError( "Couldn't find method", nsme);
    } catch (NoSuchFieldException nsfe) {
      throw new BioError( "Couldn't find field", nsfe);
    } catch (IllegalSymbolException ise) {
      throw new BioError( "Couldn't find symbol", ise);
    } catch (BioException be) {
      throw new BioError( "Couldn't create indexer", be);
    }
  }

  public Class generateViterbiClass(DP dp) {
    // viterbi recursion is:
    // score[state_i] =
    //   emission * max_j(prev_score[state_j] * transition[state_j, state_i])
    // however, in log space (as we use), this becomes:
    // score[state_i] =
    //   emission + max_j(
    //     prev_score[state_j] + transition[state_j, state_i]
    //   ))

    try {
      MarkovModel model = dp.getModel();

      String name = makeName("org.biojava.bio.dp.twohead.Viterbi", model);
      if(!classLoader.hasGeneratedClass(name)) {

      CodeClass _Object = IntrospectedCodeClass.forClass(Object.class);
      CodeClass _DP = IntrospectedCodeClass.forClass(DP.class);
      CodeClass _BackPointer = IntrospectedCodeClass.forClass(BackPointer.class);
      CodeClass _ScoreType = IntrospectedCodeClass.forClass(ScoreType.class);
      CodeClass _State_A = IntrospectedCodeClass.forClass(State [].class);
      CodeClass _Cell_A_A = IntrospectedCodeClass.forClass(Cell [][].class);
      CodeClass _CellCalculator = IntrospectedCodeClass.forClass(CellCalculator.class);

      GeneratedCodeClass clazz = new GeneratedCodeClass(
        name,
        _Object,
        new CodeClass[] {_CellCalculator},
        CodeUtils.ACC_PUBLIC | CodeUtils.ACC_STRICT
      );

      State[] states = dp.getStates();
      AlphabetIndex stateIndexer = AlphabetManager.getAlphabetIndex(states);

      CodeField terminalBP = clazz.createField(
        "terminalBP",
        _BackPointer,
        CodeUtils.ACC_PROTECTED
      );

      CodeField stateF = clazz.createField(
        "states",
        _State_A,
        CodeUtils.ACC_PROTECTED
      );

      // The fields that contain transition scores as double [].
      // This is built so that if two states shair the same transition
      // distribution, they refer to the same transition field. This gives
      // good optimizers a fighting chance to do some hard-core optimizations.
      CodeField[] transitionFields = new CodeField[states.length];
      AlphabetIndex[] transitionIndexers = new AlphabetIndex[states.length];

      // The constructor must load in the transition probabilities.
      // It uses the indexToFieldIndex to ensure that parameters are loaded
      // just once.
      GeneratedCodeMethod __init = clazz.createMethod(
        "<init>",
        CodeUtils.TYPE_VOID,
        new CodeClass[] {_DP, _ScoreType, _BackPointer },
        new String[] { "dp", "scoreType", "backPointer" },
        CodeUtils.ACC_PUBLIC
      );

      clazz.setCodeGenerator(__init, createInit(
        model, states,
        clazz, __init,
        transitionIndexers, transitionFields,
        stateF, terminalBP
      ));

      GeneratedCodeMethod initialize = clazz.createMethod(
        "initialize",
        CodeUtils.TYPE_VOID,
        new CodeClass[] {_Cell_A_A},
        new String[] { "cells" },
        CodeUtils.ACC_PUBLIC
      );
      clazz.setCodeGenerator(
        initialize,
        createVRecursion(
          true,
          model, states, stateIndexer,
          stateF, transitionFields, transitionIndexers,
          terminalBP,
          initialize
        )
      );

      GeneratedCodeMethod calcCell = clazz.createMethod(
        "calcCell",
        CodeUtils.TYPE_VOID,
        new CodeClass[] {_Cell_A_A },
        new String[] { "cells" },
        CodeUtils.ACC_PUBLIC
      );
      clazz.setCodeGenerator(
        calcCell,
        createVRecursion(
          false,
          model, states, stateIndexer,
          stateF, transitionFields, transitionIndexers,
          terminalBP,
          calcCell
        )
      );

      if(dumpToDisk == true) {
        try {
          StringBuffer fName = new StringBuffer(name);
          for(int i = 0; i < fName.length(); i++) {
            if(fName.charAt(i) == '.') {
              fName.setCharAt(i, '/');
            }
          }
          fName.append(".class");

          File dumpFile = new File(fName.substring(0));
          System.out.println("Dumping file: " + fName);
          dumpFile.getParentFile().mkdirs();

          OutputStream out = new FileOutputStream(dumpFile);
          clazz.createCode(out);
          out.close();
        } catch (FileNotFoundException fnfe) {
          throw new BioError( "Couldn't dump dp class", fnfe);
        } catch (IOException ioe) {
          throw new BioError( "Couldn't dump dp class", ioe);
        }
      }

      classLoader.defineClass(clazz);
    }
     try {
          return classLoader.loadClass(name);
        } catch (Exception e) {
          throw new BioError( "Can't find previously generated class for " + name, e);
        }
    } catch (CodeException ce) {
      throw new BioError( "Couldn't generate class", ce);
    } catch (NoSuchMethodException nsme) {
      throw new BioError( "Couldn't find method", nsme);
    } catch (NoSuchFieldException nsfe) {
      throw new BioError( "Couldn't find field", nsfe);
    } catch (IllegalSymbolException ise) {
      throw new BioError( "Couldn't find symbol", ise);
    } catch (BioException be) {
      throw new BioError( "Couldn't create indexer", be);
    }
  }

  CodeGenerator createVRecursion(
    boolean isInit,
    MarkovModel model,
    State[] states,
    AlphabetIndex stateIndex,
    CodeField stateF,
    CodeField [] transitionFields,
    AlphabetIndex[] transitionIndexers,
    CodeField terminalBP,
    GeneratedCodeMethod method
  ) throws
    NoSuchMethodException,
    NoSuchFieldException,
    IllegalSymbolException,
    CodeException
  {
    CodeClass _Cell = IntrospectedCodeClass.forClass(Cell.class);
    CodeField _Cell_score = _Cell.getFieldByName("scores");
    CodeField _Cell_backpointer = _Cell.getFieldByName("backPointers");
    CodeField _Cell_emissions = _Cell.getFieldByName("emissions");
    CodeClass _State = IntrospectedCodeClass.forClass(State.class);
    CodeClass _BackPointer = IntrospectedCodeClass.forClass(BackPointer.class);
    CodeClass _BackPointer_A = IntrospectedCodeClass.forClass(BackPointer [].class);
    CodeMethod _BackPointer_init = _BackPointer.getConstructor(new CodeClass [] {_State, _BackPointer, CodeUtils.TYPE_DOUBLE});
    CodeClass _double_A = IntrospectedCodeClass.forClass(double [].class);

    InstructionVector ccV = new InstructionVector();

    // if(isInit && state instanceof emission) {
    //   cell[0][0] = (state == Magical) ? 0.0 : NaN;
    // } else {
    //   cell[0][0].score[i] = cell[0][0].emissions[i] +
    //                         max_j(cell[adv_i_0][adv_i_1] + t_j[i])
    // }

    // cell_00 = cell[0][0];
    // cell_01 = cell[0][1];
    // cell_10 = cell[1][0];
    // cell_11 = cell[1][1];
    LocalVariable[][] cell = new LocalVariable[2][2];
    cell[0][0] = new LocalVariable(_Cell, "cell_00");
    cell[0][1] = new LocalVariable(_Cell, "cell_01");
    cell[1][0] = new LocalVariable(_Cell, "cell_10");
    cell[1][1] = new LocalVariable(_Cell, "cell_11");

    ccV.add( ByteCode.make_aload  (method.getVariable("cells")));
    ccV.add( ByteCode.make_dup    ());
    ccV.add( ByteCode.make_iconst (0));
    ccV.add( ByteCode.make_aaload ());
    ccV.add( ByteCode.make_dup    ());
    ccV.add( ByteCode.make_iconst (0));
    ccV.add( ByteCode.make_aaload ());
    ccV.add( ByteCode.make_astore (cell[0][0]));
    ccV.add( ByteCode.make_iconst (1));
    ccV.add( ByteCode.make_aaload ());
    ccV.add( ByteCode.make_astore (cell[0][1]));
    ccV.add( ByteCode.make_iconst (1));
    ccV.add( ByteCode.make_aaload ());
    ccV.add( ByteCode.make_dup    ());
    ccV.add( ByteCode.make_iconst (0));
    ccV.add( ByteCode.make_aaload ());
    ccV.add( ByteCode.make_astore (cell[1][0]));
    ccV.add( ByteCode.make_iconst (1));
    ccV.add( ByteCode.make_aaload ());
    ccV.add( ByteCode.make_astore (cell[1][1]));

    // score_00 = cell[0][0].score;
    // score_01 = cell[0][1].score;
    // score_10 = cell[1][0].score;
    // score_11 = cell[1][1].score;

    LocalVariable[][] score = new LocalVariable[2][2];
    score[0][0] = new LocalVariable(_double_A, "score_00");
    score[0][1] = new LocalVariable(_double_A, "score_01");
    score[1][0] = new LocalVariable(_double_A, "score_10");
    score[1][1] = new LocalVariable(_double_A, "score_11");
    ccV.add( ByteCode.make_aload    (cell[0][0]));
    ccV.add( ByteCode.make_getfield ( _Cell_score));
    ccV.add( ByteCode.make_astore   (score[0][0]));
    ccV.add( ByteCode.make_aload    (cell[0][1]));
    ccV.add( ByteCode.make_getfield ( _Cell_score));
    ccV.add( ByteCode.make_astore   (score[0][1]));
    ccV.add( ByteCode.make_aload    (cell[1][0]));
    ccV.add( ByteCode.make_getfield ( _Cell_score));
    ccV.add( ByteCode.make_astore   (score[1][0]));
    ccV.add( ByteCode.make_aload    (cell[1][1]));
    ccV.add( ByteCode.make_getfield ( _Cell_score));
    ccV.add( ByteCode.make_astore   (score[1][1]));

    // backpointer_00 = cell[0][0].backpointer;
    // backpointer_01 = cell[0][1].backpointer;
    // backpointer_10 = cell[1][0].backpointer;
    // backpointer_11 = cell[1][1].backpointer;

    LocalVariable[][] backpointer = new LocalVariable[2][2];
    backpointer[0][0] = new LocalVariable(_BackPointer_A, "backpointer_00");
    backpointer[0][1] = new LocalVariable(_BackPointer_A, "backpointer_01");
    backpointer[1][0] = new LocalVariable(_BackPointer_A, "backpointer_10");
    backpointer[1][1] = new LocalVariable(_BackPointer_A, "backpointer_11");
    ccV.add( ByteCode.make_aload    (cell[0][0]));
    ccV.add( ByteCode.make_getfield ( _Cell_backpointer));
    ccV.add( ByteCode.make_astore   (backpointer[0][0]));
    ccV.add( ByteCode.make_aload    (cell[0][1]));
    ccV.add( ByteCode.make_getfield ( _Cell_backpointer));
    ccV.add( ByteCode.make_astore   (backpointer[0][1]));
    ccV.add( ByteCode.make_aload    (cell[1][0]));
    ccV.add( ByteCode.make_getfield ( _Cell_backpointer));
    ccV.add( ByteCode.make_astore   (backpointer[1][0]));
    ccV.add( ByteCode.make_aload    (cell[1][1]));
    ccV.add( ByteCode.make_getfield ( _Cell_backpointer));
    ccV.add( ByteCode.make_astore   (backpointer[1][1]));

    LocalVariable emissions = new LocalVariable(_double_A, "emissions");
    ccV.add( ByteCode.make_aload    (cell[0][0] ));
    ccV.add( ByteCode.make_getfield (_Cell_emissions));
    ccV.add( ByteCode.make_astore   (emissions));

    LocalVariable max = new LocalVariable(CodeUtils.TYPE_DOUBLE, "max");
    LocalVariable max_j = new LocalVariable(CodeUtils.TYPE_INT, "max_j");
    for(int i = 0; i < states.length; i++) {
      State state = states[i];
      InstructionVector stateV = new InstructionVector();
      if(isInit && state instanceof EmissionState) {
        stateV.add( ByteCode.make_aload   (score[0][0]));
        stateV.add( ByteCode.make_iconst  (i));
        if(state instanceof MagicalState) {
          stateV.add( ByteCode.make_dconst   (0.0));
          stateV.add( ByteCode.make_aload    (backpointer[0][0]));
          stateV.add( ByteCode.make_iconst   (i));
          stateV.add( ByteCode.make_aload    (method.getThis()));
          stateV.add( ByteCode.make_getfield (terminalBP));
          stateV.add( ByteCode.make_aastore  ());
        } else {
          stateV.add( ByteCode.make_dconst   (Double.NaN));
        }
        stateV.add( ByteCode.make_dastore ());
      } else {
        int[] advance = getAdvance(state);
        FiniteAlphabet trans = model.transitionsFrom(state);

        // find max/argmax of t_j[i] + v[j]
        // make a max pipeline
        Iterator each_j = trans.iterator();
        State state_j;
        int state_jIndx;

        // first state primes the pump
        stateV.add( ByteCode.make_dconst (Double.NEGATIVE_INFINITY));
        stateV.add( ByteCode.make_dstore (max));
        stateV.add( ByteCode.make_iconst (-1));
        stateV.add( ByteCode.make_istore (max_j));

        // then process each state
        while(each_j.hasNext()) {
          state_j = (State) each_j.next();
          state_jIndx = stateIndex.indexForSymbol(state_j);
          stateV.add( createTransitionLastSum(
            method,
            transitionIndexers[state_jIndx].indexForSymbol(state),
            state_jIndx,
            transitionFields,
            getAdvanced(score, advance)
          ));
          stateV.add( ByteCode.make_dup2  ());
          stateV.add( ByteCode.make_dload (max));
          stateV.add( ByteCode.make_dcmpl ());

          // if dcmpl is 1, we need to store the new max & imax. If either
          // is NaN, keep the current max.
          InstructionVector saveNewMax = new InstructionVector();
          saveNewMax.add( ByteCode.make_dstore (max));
          saveNewMax.add( ByteCode.make_iconst (state_jIndx));
          saveNewMax.add( ByteCode.make_istore (max_j));

          // if they are equal or max is greater or either is NaN
          // dump current value
          InstructionVector useOldMax = new InstructionVector();
          useOldMax.add( ByteCode.make_pop2());

          stateV.add( new IfExpression(
            ByteCode.op_ifge, // branch if int on stack is >= 0
            saveNewMax,
            useOldMax
          ));
        }

        // if (max_j == -1)
        //   score[i] = NaN
        //   bp[i] = null
        // else
        //   sum = emissions[i] + max
        //   score[i] = sum
        //   bp[i] = new BackPointer(state[i], bp[adv_0][adv_1], sum)
        // endif

        // (max == max) && (max_j != -1)
        stateV.add( ByteCode.make_iload  (max_j));

        InstructionVector ifNoIn = new InstructionVector();
        InstructionVector ifGotIn = new InstructionVector();

        ifNoIn.add( ByteCode.make_aload       (score[0][0]));
        ifNoIn.add( ByteCode.make_iconst      (i));
        ifNoIn.add( ByteCode.make_dconst      (Double.NaN));
        ifNoIn.add( ByteCode.make_dastore     ());
        ifNoIn.add( ByteCode.make_aload       (backpointer[0][0]));
        ifNoIn.add( ByteCode.make_iconst      (i));
        ifNoIn.add( ByteCode.make_aconst_null ());
        ifNoIn.add( ByteCode.make_aastore     ());

        // score[i] = emissions[i] + max if emitting state, else score[i] = max
        ifGotIn.add( ByteCode.make_aload   (score[0][0]));
        ifGotIn.add( ByteCode.make_iconst  (i));
        ifGotIn.add( ByteCode.make_dload   (max));
        if(state instanceof EmissionState) {
          ifGotIn.add( ByteCode.make_aload   (emissions));
          ifGotIn.add( ByteCode.make_iconst  (i));
          ifGotIn.add( ByteCode.make_daload  ());
          ifGotIn.add( ByteCode.make_dadd    ());
        }
        ifGotIn.add( ByteCode.make_dastore ());

        // backpointer[i] = new BackPointer(
        //  state[i],
        //  backPointer[adv_0][adv_1],
        //  score
        // );

        // backpointer[i] =
        ifGotIn.add( ByteCode.make_aload    (backpointer[0][0]));
        ifGotIn.add( ByteCode.make_iconst   (i));

        // new BackPointer - dup for <init> invoke
        ifGotIn.add( ByteCode.make_new (_BackPointer));
        ifGotIn.add( ByteCode.make_dup ());

        // state[i]
        ifGotIn.add( ByteCode.make_aload    (method.getThis()));
        ifGotIn.add( ByteCode.make_getfield (stateF));
        ifGotIn.add( ByteCode.make_iconst   (i));
        ifGotIn.add( ByteCode.make_aaload   ());

        // backpointer[adv_0][adv_1] [max_j]
        ifGotIn.add( ByteCode.make_aload  (getAdvanced(backpointer, advance)));
        ifGotIn.add( ByteCode.make_iload  (max_j));
        ifGotIn.add( ByteCode.make_aaload ());

        // score[i]
        ifGotIn.add( ByteCode.make_aload  (score[0][0]));
        ifGotIn.add( ByteCode.make_iconst (i));
        ifGotIn.add( ByteCode.make_daload ());

        // backpointer.<init>
        ifGotIn.add( ByteCode.make_invokespecial( _BackPointer_init));

        // store backpointer
        ifGotIn.add( ByteCode.make_aastore ());

        stateV.add( new IfExpression(
          ByteCode.op_ifge, // ifge means max_j >= 0
          ifGotIn,
          ifNoIn
        ));
      }
     ccV.add(stateV);
    }

    // dump out scores & backpointers
    /*LocalVariable sc = new LocalVariable(CodeUtils.TYPE_DOUBLE, "score_i");
    LocalVariable bp = new LocalVariable(_BackPointer, "bp_i");
    for(int i = 0; i < states.length; i++) {
      ccV.add( ByteCode.make_aload  (score[0][0]));
      ccV.add( ByteCode.make_iconst (i));
      ccV.add( ByteCode.make_daload ());
      ccV.add( ByteCode.make_dstore (sc));

      ccV.add( ByteCode.make_aload  (backpointer[0][0]));
      ccV.add( ByteCode.make_iconst (i));
      ccV.add( ByteCode.make_aaload ());
      ccV.add( ByteCode.make_astore (bp));
    }
    */
    ccV.add( ByteCode.make_return ());

    return ccV;
  }

  private CodeGenerator createFRecursion(
    boolean isInit,
    MarkovModel model,
    State[] states,
    AlphabetIndex stateIndex,
    CodeField stateF,
    CodeField [] transitionFields,
    AlphabetIndex[] transitionIndexers,
    GeneratedCodeMethod method
  ) throws
    NoSuchMethodException,
    NoSuchFieldException,
    IllegalSymbolException,
    CodeException
  {
    CodeClass _Cell = IntrospectedCodeClass.forClass(Cell.class);
    CodeField _Cell_score = _Cell.getFieldByName("scores");
    CodeField _Cell_emissions = _Cell.getFieldByName("emissions");
    CodeClass _double_A = IntrospectedCodeClass.forClass(double [].class);
    CodeClass _Math = IntrospectedCodeClass.forClass(Math.class);
    CodeMethod _Math_exp = _Math.getMethod("exp", new CodeClass[] { CodeUtils.TYPE_DOUBLE });
    CodeMethod _Math_log = _Math.getMethod("log", new CodeClass[] { CodeUtils.TYPE_DOUBLE });

    InstructionVector ccV = new InstructionVector();

    ccV.add( debug(message("Retrieving local variables")));

    // if(isInit && state instanceof emission) {
    //   cell[0][0] = (state == Magical) ? 0.0 : NaN;
    // } else {
    //   cell[0][0].score[i] = cell[0][0].emissions[i] +
    //                         sum_j(cell[adv_i_0][adv_i_1] + t_j[i])
    // }

    // cell_00 = cell[0][0];
    // cell_01 = cell[0][1];
    // cell_10 = cell[1][0];
    // cell_11 = cell[1][1];
    LocalVariable[][] cell = new LocalVariable[2][2];
    cell[0][0] = new LocalVariable(_Cell, "cell_00");
    cell[0][1] = new LocalVariable(_Cell, "cell_01");
    cell[1][0] = new LocalVariable(_Cell, "cell_10");
    cell[1][1] = new LocalVariable(_Cell, "cell_11");

    ccV.add( ByteCode.make_aload  (method.getVariable("cells")));
    ccV.add( ByteCode.make_dup    ());
    ccV.add( ByteCode.make_iconst (0));
    ccV.add( ByteCode.make_aaload ());
    ccV.add( ByteCode.make_dup    ());
    ccV.add( ByteCode.make_iconst (0));
    ccV.add( ByteCode.make_aaload ());
    ccV.add( ByteCode.make_astore (cell[0][0]));
    ccV.add( ByteCode.make_iconst (1));
    ccV.add( ByteCode.make_aaload ());
    ccV.add( ByteCode.make_astore (cell[0][1]));
    ccV.add( ByteCode.make_iconst (1));
    ccV.add( ByteCode.make_aaload ());
    ccV.add( ByteCode.make_dup    ());
    ccV.add( ByteCode.make_iconst (0));
    ccV.add( ByteCode.make_aaload ());
    ccV.add( ByteCode.make_astore (cell[1][0]));
    ccV.add( ByteCode.make_iconst (1));
    ccV.add( ByteCode.make_aaload ());
    ccV.add( ByteCode.make_astore (cell[1][1]));

    // score_00 = cell[0][0].score;
    // score_01 = cell[0][1].score;
    // score_10 = cell[1][0].score;
    // score_11 = cell[1][1].score;

    LocalVariable[][] score = new LocalVariable[2][2];
    score[0][0] = new LocalVariable(_double_A, "score_00");
    score[0][1] = new LocalVariable(_double_A, "score_01");
    score[1][0] = new LocalVariable(_double_A, "score_10");
    score[1][1] = new LocalVariable(_double_A, "score_11");
    ccV.add( ByteCode.make_aload    (cell[0][0]));
    ccV.add( ByteCode.make_getfield ( _Cell_score));
    ccV.add( ByteCode.make_astore   (score[0][0]));
    ccV.add( ByteCode.make_aload    (cell[0][1]));
    ccV.add( ByteCode.make_getfield ( _Cell_score));
    ccV.add( ByteCode.make_astore   (score[0][1]));
    ccV.add( ByteCode.make_aload    (cell[1][0]));
    ccV.add( ByteCode.make_getfield ( _Cell_score));
    ccV.add( ByteCode.make_astore   (score[1][0]));
    ccV.add( ByteCode.make_aload    (cell[1][1]));
    ccV.add( ByteCode.make_getfield ( _Cell_score));
    ccV.add( ByteCode.make_astore   (score[1][1]));

    LocalVariable emissions = new LocalVariable(_double_A, "emissions");
    ccV.add( ByteCode.make_aload    (cell[0][0] ));
    ccV.add( ByteCode.make_getfield (_Cell_emissions));
    ccV.add( ByteCode.make_astore   (emissions));

    LocalVariable max = new LocalVariable(CodeUtils.TYPE_DOUBLE, "max");
    for(int i = 0; i < states.length; i++) {
      State state = states[i];
      InstructionVector stateV = new InstructionVector();
      // we need to push score & i onto the stack so that after finding the sum
      // we can just push it back into the array
      stateV.add( ByteCode.make_aload   (score[0][0]));
      stateV.add( ByteCode.make_iconst  (i));
      if(isInit && state instanceof EmissionState) {
        if(state instanceof MagicalState) {
          stateV.add( ByteCode.make_dconst   (0.0));
        } else {
          stateV.add( ByteCode.make_dconst   (Double.NaN));
        }
      } else {
        int[] advance = getAdvance(state);
        FiniteAlphabet trans = model.transitionsFrom(state);

        LocalVariable j_scores = getAdvanced(score, advance);

        if(trans.size() == 1) {
          State state_j = (State) trans.iterator().next();
          int state_jIndx = stateIndex.indexForSymbol(state_j);

          // only one source. F[i] = trans_j[i] + F[j] + e[i]
          stateV.add( ByteCode.make_aload    (method.getThis()));
          stateV.add( ByteCode.make_getfield (transitionFields[state_jIndx]));
          stateV.add( ByteCode.make_iconst   (transitionIndexers[state_jIndx].indexForSymbol(state)));
          stateV.add( ByteCode.make_daload   ());
          stateV.add( ByteCode.make_aload    (j_scores));
          stateV.add( ByteCode.make_iconst   (state_jIndx));
          stateV.add( ByteCode.make_daload   ());
          stateV.add( ByteCode.make_dadd     ());
          if(state instanceof EmissionState) {
            stateV.add( ByteCode.make_aload  (emissions));
            stateV.add( ByteCode.make_iconst (i));
            stateV.add( ByteCode.make_daload ());
            stateV.add( ByteCode.make_dadd   ());
          }
        } else {
          // f[i] = emission[i] + log( sum_j [
          //   exp( (j_scores[j] - max_j_score) + t_j[i] + e)
          // ]) + max_j_score

          // find max j_scores
          stateV.add( ByteCode.make_dconst (Double.NEGATIVE_INFINITY));
          stateV.add( ByteCode.make_dstore (max));

          Iterator each_j = trans.iterator();
          while(each_j.hasNext()) {
            State state_j = (State) each_j.next();
            int state_jIndx = stateIndex.indexForSymbol(state_j);
            stateV.add( ByteCode.make_aload    (j_scores));
            stateV.add( ByteCode.make_iconst   (state_jIndx));
            stateV.add( ByteCode.make_daload   ());
            stateV.add( ByteCode.make_dup2     ());
            stateV.add( ByteCode.make_dload    (max));
            stateV.add( ByteCode.make_dcmpl    ()); // puts -1 if max > this one

            InstructionVector ifLargerThanMax = new InstructionVector();
            ifLargerThanMax.add( ByteCode.make_dstore (max));

            InstructionVector ifSmallerThanMax = new InstructionVector();
            ifSmallerThanMax.add( ByteCode.make_pop2());

            stateV.add( new IfExpression(
              ByteCode.op_ifge, // branch if int on stack >= 0
              ifLargerThanMax,
              ifSmallerThanMax
            ));
          }

          LocalVariable maxDB = new LocalVariable(CodeUtils.TYPE_DOUBLE, "max");
          InstructionVector dbi = new InstructionVector();
          dbi.add( ByteCode.make_dload (max));
          dbi.add( ByteCode.make_dstore (maxDB));
          dbi.add( message("Max " + i + " = ", maxDB));
          stateV.add( debug(dbi));

          // log sum of exponents - prime the sum with zero
          stateV.add( ByteCode.make_dconst (0.0));

          each_j = trans.iterator();
          while(each_j.hasNext()) {
            State state_j = (State) each_j.next();
            int state_jIndx = stateIndex.indexForSymbol(state_j);

            // score
            stateV.add( ByteCode.make_aload  (j_scores));
            stateV.add( ByteCode.make_iconst (state_jIndx));
            stateV.add( ByteCode.make_daload ());

            // is score NaN?
            stateV.add( ByteCode.make_dup2());
            stateV.add( ByteCode.make_dup2());
            stateV.add( ByteCode.make_dcmpl());

            InstructionVector scoreNotNaN = new InstructionVector();
            // load max & subtract it from score
            // (j_score - max)
            scoreNotNaN.add( ByteCode.make_dload  (max));
            scoreNotNaN.add( ByteCode.make_dsub   ());

            // exp( (j_score - max) + transition)
            scoreNotNaN.add( ByteCode.make_aload        (method.getThis()));
            scoreNotNaN.add( ByteCode.make_getfield     (transitionFields[state_jIndx]));
            scoreNotNaN.add( ByteCode.make_iconst       (transitionIndexers[state_jIndx].indexForSymbol(state)));
            scoreNotNaN.add( ByteCode.make_daload       ());
            scoreNotNaN.add( ByteCode.make_dadd         ());
            scoreNotNaN.add( ByteCode.make_invokestatic (_Math_exp));

            // sum this and current sum
            scoreNotNaN.add( ByteCode.make_dadd());

            InstructionVector scoreIsNaN = new InstructionVector();
            scoreIsNaN.add( ByteCode.make_pop2());

            stateV.add( new IfExpression(
              ByteCode.op_ifge,
              scoreNotNaN,
              scoreIsNaN
            ));
          }

          // log sum
          stateV.add( ByteCode.make_invokestatic (_Math_log));

          if(state instanceof EmissionState) {
            // emissions[i]
            stateV.add( ByteCode.make_aload    (emissions));
            stateV.add( ByteCode.make_iconst   (i));
            stateV.add( ByteCode.make_daload   ());

            // sum emission with j sum
            stateV.add( ByteCode.make_dadd     ());
          }

          // lastly add on a factor of max - added here for maximum numerical stability
          stateV.add( ByteCode.make_dload (max));
          stateV.add( ByteCode.make_dadd  ());
        }
      }
      // store score on stack into scores array
      stateV.add( ByteCode.make_dastore ());
      ccV.add( stateV );
    }

    // dump out state scores
    LocalVariable sc = new LocalVariable(CodeUtils.TYPE_DOUBLE, "score_i");
    InstructionVector dbi = new InstructionVector();
    for(int i = 0; i < states.length; i++) {
      dbi.add( ByteCode.make_aload  (score[0][0]));
      dbi.add( ByteCode.make_iconst (i));
      dbi.add( ByteCode.make_daload ());
      dbi.add( ByteCode.make_dstore (sc));
      dbi.add( message("Score " + i + " = ", sc));
    }
    ccV.add( debug(dbi));

    ccV.add( ByteCode.make_return ());

    return ccV;
  }

  private CodeGenerator createBRecursion(
    boolean isInit,
    MarkovModel model,
    State[] states,
    AlphabetIndex stateIndex,
    CodeField stateF,
    CodeField [] transitionFields,
    AlphabetIndex[] transitionIndexers,
    GeneratedCodeMethod method
  ) throws
    NoSuchMethodException,
    NoSuchFieldException,
    IllegalSymbolException,
    CodeException
  {
    CodeClass _Cell = IntrospectedCodeClass.forClass(Cell.class);
    CodeField _Cell_score = _Cell.getFieldByName("scores");
    CodeField _Cell_emissions = _Cell.getFieldByName("emissions");
    CodeClass _double_A = IntrospectedCodeClass.forClass(double [].class);
    CodeClass _Math = IntrospectedCodeClass.forClass(Math.class);
    CodeMethod _Math_exp = _Math.getMethod("exp", new CodeClass[] { CodeUtils.TYPE_DOUBLE });
    CodeMethod _Math_log = _Math.getMethod("log", new CodeClass[] { CodeUtils.TYPE_DOUBLE });

    InstructionVector ccV = new InstructionVector();

    ccV.add( debug(message("Retrieving local variables")));

    // if(isInit && state instanceof emission) {
    //   cell[0][0] = (state == Magical) ? 0.0 : NaN;
    // } else {
    //   cell[0][0].score[i] = cell[0][0].emissions[i] +
    //                         sum_j(cell[adv_i_0][adv_i_1] + t_j[i])
    // }

    // cell_00 = cell[0][0];
    // cell_01 = cell[0][1];
    // cell_10 = cell[1][0];
    // cell_11 = cell[1][1];
    LocalVariable[][] cell = new LocalVariable[2][2];
    cell[0][0] = new LocalVariable(_Cell, "cell_00");
    cell[0][1] = new LocalVariable(_Cell, "cell_01");
    cell[1][0] = new LocalVariable(_Cell, "cell_10");
    cell[1][1] = new LocalVariable(_Cell, "cell_11");

    ccV.add( ByteCode.make_aload  (method.getVariable("cells")));
    ccV.add( ByteCode.make_dup    ());
    ccV.add( ByteCode.make_iconst (0));
    ccV.add( ByteCode.make_aaload ());
    ccV.add( ByteCode.make_dup    ());
    ccV.add( ByteCode.make_iconst (0));
    ccV.add( ByteCode.make_aaload ());
    ccV.add( ByteCode.make_astore (cell[0][0]));
    ccV.add( ByteCode.make_iconst (1));
    ccV.add( ByteCode.make_aaload ());
    ccV.add( ByteCode.make_astore (cell[0][1]));
    ccV.add( ByteCode.make_iconst (1));
    ccV.add( ByteCode.make_aaload ());
    ccV.add( ByteCode.make_dup    ());
    ccV.add( ByteCode.make_iconst (0));
    ccV.add( ByteCode.make_aaload ());
    ccV.add( ByteCode.make_astore (cell[1][0]));
    ccV.add( ByteCode.make_iconst (1));
    ccV.add( ByteCode.make_aaload ());
    ccV.add( ByteCode.make_astore (cell[1][1]));

    // score_00 = cell[0][0].score;
    // score_01 = cell[0][1].score;
    // score_10 = cell[1][0].score;
    // score_11 = cell[1][1].score;

    LocalVariable[][] score = new LocalVariable[2][2];
    score[0][0] = new LocalVariable(_double_A, "score_00");
    score[0][1] = new LocalVariable(_double_A, "score_01");
    score[1][0] = new LocalVariable(_double_A, "score_10");
    score[1][1] = new LocalVariable(_double_A, "score_11");
    ccV.add( ByteCode.make_aload    (cell[0][0]));
    ccV.add( ByteCode.make_getfield ( _Cell_score));
    ccV.add( ByteCode.make_astore   (score[0][0]));
    ccV.add( ByteCode.make_aload    (cell[0][1]));
    ccV.add( ByteCode.make_getfield ( _Cell_score));
    ccV.add( ByteCode.make_astore   (score[0][1]));
    ccV.add( ByteCode.make_aload    (cell[1][0]));
    ccV.add( ByteCode.make_getfield ( _Cell_score));
    ccV.add( ByteCode.make_astore   (score[1][0]));
    ccV.add( ByteCode.make_aload    (cell[1][1]));
    ccV.add( ByteCode.make_getfield ( _Cell_score));
    ccV.add( ByteCode.make_astore   (score[1][1]));

    LocalVariable[][] emission = new LocalVariable[2][2];
    emission[0][1] = new LocalVariable(_double_A, "emission_01");
    emission[1][0] = new LocalVariable(_double_A, "emission_10");
    emission[1][1] = new LocalVariable(_double_A, "emission_11");
    ccV.add( ByteCode.make_aload    (cell[0][1] ));
    ccV.add( ByteCode.make_getfield (_Cell_emissions));
    ccV.add( ByteCode.make_astore   (emission[0][1]));
    ccV.add( ByteCode.make_aload    (cell[1][0] ));
    ccV.add( ByteCode.make_getfield (_Cell_emissions));
    ccV.add( ByteCode.make_astore   (emission[1][0]));
    ccV.add( ByteCode.make_aload    (cell[1][1] ));
    ccV.add( ByteCode.make_getfield (_Cell_emissions));
    ccV.add( ByteCode.make_astore   (emission[1][1]));

    LocalVariable max = new LocalVariable(CodeUtils.TYPE_DOUBLE, "max");
    for(int i = states.length-1; i >= 0 ; i--) {
      State state = states[i];
      ccV.add( debug(message("Calculating for state " + i + " " + state.getName())));
      InstructionVector stateV = new InstructionVector();
      // we need to push score & i onto the stack so that after finding the sum
      // we can just push it back into the array
      stateV.add( ByteCode.make_aload   (score[0][0]));
      stateV.add( ByteCode.make_iconst  (i));
      if(isInit && state instanceof EmissionState) {
        stateV.add( debug(message("initalizing")));
        if(state instanceof MagicalState) {
          stateV.add( debug(message("magical")));
          stateV.add( ByteCode.make_dconst   (0.0));
        } else {
          stateV.add( debug(message("mundane")));
          stateV.add( ByteCode.make_dconst   (Double.NaN));
        }
      } else {
        FiniteAlphabet trans = model.transitionsFrom(state);

        if(trans.size() == 1) {
          stateV.add( debug(message("single-source optimization")));

          State state_j = (State) trans.iterator().next();
          int state_jIndx = stateIndex.indexForSymbol(state_j);

          // work out advance
          // only one source. B[i] = trans_j[i] + B[j] + e[j]
          // only add e[j] if emitting state
          int [] advance = getAdvance(state_j);

          stateV.add( ByteCode.make_aload    (method.getThis()));
          stateV.add( ByteCode.make_getfield (transitionFields[i]));
          stateV.add( ByteCode.make_iconst   (transitionIndexers[i].indexForSymbol(state_j)));
          stateV.add( ByteCode.make_daload   ());
          stateV.add( ByteCode.make_aload    (getAdvanced(score, advance)));
          stateV.add( ByteCode.make_iconst   (state_jIndx));
          stateV.add( ByteCode.make_daload   ());
          stateV.add( ByteCode.make_dadd     ());
          if(state_j instanceof EmissionState) {
            stateV.add( ByteCode.make_aload  (getAdvanced(emission, advance)));
            stateV.add( ByteCode.make_iconst (i));
            stateV.add( ByteCode.make_daload ());
            stateV.add( ByteCode.make_dadd   ());
          }

          LocalVariable tmp = new LocalVariable(CodeUtils.TYPE_DOUBLE, "tmp");
          InstructionVector dbi = new InstructionVector();
          dbi.add( ByteCode.make_dup2());
          dbi.add( ByteCode.make_dstore(tmp));
          dbi.add( message("got sum of ", tmp));
          stateV.add( debug(dbi));

        } else {
          // B[i] = log( sum_j [
          //   exp( (j_scores[j] - max_j_score) + t_j[i] + e[j])
          // ]) + max_j_score

          stateV.add( debug(message("full recursion")));

          // find max j_scores
          stateV.add( ByteCode.make_dconst (Double.NEGATIVE_INFINITY));
          stateV.add( ByteCode.make_dstore (max));

          Iterator each_j = trans.iterator();
          while(each_j.hasNext()) {
            State state_j = (State) each_j.next();
            int state_jIndx = stateIndex.indexForSymbol(state_j);
            int[] advance = getAdvance(state_j);

            stateV.add( ByteCode.make_aload    (getAdvanced(score, advance)));
            stateV.add( ByteCode.make_iconst   (state_jIndx));
            stateV.add( ByteCode.make_daload   ());
            stateV.add( ByteCode.make_dup2     ());
            stateV.add( ByteCode.make_dload    (max));
            stateV.add( ByteCode.make_dcmpl    ()); // puts -1 if max > this one

            InstructionVector ifLargerThanMax = new InstructionVector();
            ifLargerThanMax.add( debug(message("Larger")));
            ifLargerThanMax.add( debug(message("  max: ", max)));
            ifLargerThanMax.add( ByteCode.make_dstore (max));

            InstructionVector ifSmallerThanMax = new InstructionVector();
            ifSmallerThanMax.add( debug(message("Smaller")));
            ifSmallerThanMax.add( debug(message("  max: ", max)));
            ifSmallerThanMax.add( ByteCode.make_pop2 ());

            stateV.add( new IfExpression(
              ByteCode.op_ifge, // branch if int on stack >= 0
              ifLargerThanMax,
              ifSmallerThanMax
            ));
          }

          stateV.add( debug(message("Taking out factor ", max)));

          // sum logs of exponents - prime sum with zero
          stateV.add( ByteCode.make_dconst (0.0));

          each_j = trans.iterator();
          while(each_j.hasNext()) {
            State state_j = (State) each_j.next();
            int state_jIndx = stateIndex.indexForSymbol(state_j);
            int[] advance = getAdvance(state_j);
            // score
            stateV.add( ByteCode.make_aload  (getAdvanced(score, advance)));
            stateV.add( ByteCode.make_iconst (state_jIndx));
            stateV.add( ByteCode.make_daload ());

            // is score NaN?
            stateV.add( ByteCode.make_dup2());
            stateV.add( ByteCode.make_dup2());
            stateV.add( ByteCode.make_dcmpl());
            // stack now looks like sum, score_j, score_j isNaN
            // - boolean popped off during conditional branch leaving
            //   sum, score_j

            InstructionVector scoreNotNaN = new InstructionVector();
            // load max & subtract it from score
            scoreNotNaN.add( ByteCode.make_dload  (max));
            scoreNotNaN.add( ByteCode.make_dsub   ());
            // sum, score_j - max

            // (j_score - max) + transition + emission[j]) // only add emission if emiting state
            scoreNotNaN.add( ByteCode.make_aload        (method.getThis()));
            scoreNotNaN.add( ByteCode.make_getfield     (transitionFields[i]));
            scoreNotNaN.add( ByteCode.make_iconst       (transitionIndexers[i].indexForSymbol(state_j)));
            scoreNotNaN.add( ByteCode.make_daload       ());
            scoreNotNaN.add( ByteCode.make_dadd         ());
            if(state_j instanceof EmissionState) {
              scoreNotNaN.add( ByteCode.make_aload      (getAdvanced(score, advance)));
              scoreNotNaN.add( ByteCode.make_iconst     (state_jIndx));
              scoreNotNaN.add( ByteCode.make_daload     ());
              scoreNotNaN.add( ByteCode.make_dadd       ());
            }
            scoreNotNaN.add( ByteCode.make_invokestatic (_Math_exp));

            // sum this and current sum
            scoreNotNaN.add( ByteCode.make_dadd());

            InstructionVector scoreIsNaN = new InstructionVector();
            scoreIsNaN.add( ByteCode.make_pop2());

            stateV.add( new IfExpression(
              ByteCode.op_ifge,
              scoreNotNaN,
              scoreIsNaN
            ));
          }

          // log sum
          stateV.add( ByteCode.make_invokestatic (_Math_log));

          // lastly add on a factor of max - added here for maximum numerical stability
          stateV.add( ByteCode.make_dload (max));
          stateV.add( ByteCode.make_dadd  ());
        }
      }
      // store score on stack into scores array
      stateV.add( ByteCode.make_dastore ());
      ccV.add( stateV );
    }

    // dump out state scores
    LocalVariable sc = new LocalVariable(CodeUtils.TYPE_DOUBLE, "score_i");
    InstructionVector dbi = new InstructionVector();
    for(int i = 0; i < states.length; i++) {
      dbi.add( ByteCode.make_aload  (score[0][0]));
      dbi.add( ByteCode.make_iconst (i));
      dbi.add( ByteCode.make_daload ());
      dbi.add( ByteCode.make_dstore (sc));
      dbi.add( message("Score " + i + " = ", sc));
    }
    ccV.add( debug(dbi));

    ccV.add( ByteCode.make_return ());

    return ccV;
  }

  private CodeGenerator createInit(
    MarkovModel model, State[] states,
    GeneratedCodeClass clazz,
    GeneratedCodeMethod __init,
    AlphabetIndex[] transitionIndexers,
    CodeField[] transitionFields,
    CodeField stateF,
    CodeField terminalBP
  ) throws
    NoSuchFieldException, NoSuchMethodException, CodeException,
    IllegalSymbolException, BioException
  {
      CodeClass _Object = IntrospectedCodeClass.forClass(Object.class);
      CodeMethod _Object_init = _Object.getConstructor(CodeUtils.EMPTY_LIST);
      CodeClass _Symbol = IntrospectedCodeClass.forClass(Symbol.class);
      CodeClass _State = IntrospectedCodeClass.forClass(State.class);
      CodeClass _State_A = IntrospectedCodeClass.forClass(State [].class);
      CodeClass _double_A = IntrospectedCodeClass.forClass(double [].class);
      CodeClass _DP = IntrospectedCodeClass.forClass(DP.class);
      CodeMethod _DP_getStates = _DP.getMethod("getStates", CodeUtils.EMPTY_LIST);
      CodeMethod _DP_getModel = _DP.getMethod("getModel", CodeUtils.EMPTY_LIST);
      CodeClass _Distribution = IntrospectedCodeClass.forClass(Distribution.class);
      CodeMethod _Distribution_getAlphabet = _Distribution.getMethod("getAlphabet", CodeUtils.EMPTY_LIST);
      CodeClass _ScoreType = IntrospectedCodeClass.forClass(ScoreType.class);
      CodeMethod _ScoreType_calculateScore = _ScoreType.getMethod("calculateScore", new CodeClass[] {_Distribution, _Symbol});
      CodeClass _MarkovModel = IntrospectedCodeClass.forClass(MarkovModel.class);
      CodeMethod _MarkovModel_getWeights = _MarkovModel.getMethod("getWeights", new CodeClass[] {_State});
      CodeClass _FiniteAlphabet = IntrospectedCodeClass.forClass(FiniteAlphabet.class);
      CodeMethod _FiniteAlphabet_size = _FiniteAlphabet.getMethod("size", CodeUtils.EMPTY_LIST);
      CodeClass _Math = IntrospectedCodeClass.forClass(Math.class);
      CodeMethod _Math_log = _Math.getMethod("log", new CodeClass[] { CodeUtils.TYPE_DOUBLE });

      int[] indexToFieldIndex = new int[states.length];
      Map distToIndx = new HashMap();
      for(int i = 0; i < states.length; i++) {
        State s = states[i];
        Distribution dist = model.getWeights(s);
        Integer indxI = (Integer) distToIndx.get(dist);
        if(indxI == null) {
          indxI = new Integer(i);
          distToIndx.put(dist, indxI);
          transitionFields[i] = clazz.createField(
            "t_" + i,
            _double_A,
            CodeUtils.ACC_PROTECTED
          );
          indexToFieldIndex[i] = i;
        } else {
          int indx = indxI.intValue();
          transitionFields[i] = transitionFields[indx];
          indexToFieldIndex[i] = indexToFieldIndex[indx];
        }
      }

      InstructionVector initG = new InstructionVector();
      // invoke super()
      initG.add( ByteCode.make_aload         (__init.getThis()));
      initG.add( ByteCode.make_invokespecial (_Object_init));

      // load up the dp object.
      // Store the states array, and HMM
      LocalVariable statesLV = new LocalVariable(_State_A, "states");
      LocalVariable modelLV = new LocalVariable(_MarkovModel, "model");
      initG.add( ByteCode.make_aload         (__init.getVariable("dp")));
      initG.add( ByteCode.make_dup           ());
      initG.add( ByteCode.make_invokevirtual (_DP_getStates));
      initG.add( ByteCode.make_astore        (statesLV));
      initG.add( ByteCode.make_invokevirtual (_DP_getModel));
      initG.add( ByteCode.make_astore        (modelLV));
      initG.add( ByteCode.make_aload         (__init.getThis()));
      initG.add( ByteCode.make_aload         (statesLV));
      initG.add( ByteCode.make_putfield      (stateF));

      // store the backPointer thing in terminalBP
      if(terminalBP != null) {
        initG.add( ByteCode.make_aload    (__init.getThis()));
        initG.add( ByteCode.make_aload    (__init.getVariable("backPointer")));
        initG.add( ByteCode.make_putfield (terminalBP));
      }

      LocalVariable distLV = new LocalVariable(_Distribution, "dist");

      // load in the transition probabilities to the transition fields
      for(int i = 0; i < transitionFields.length; i++) {
        // if this field reference is the first one for this distribution
        if(indexToFieldIndex[i] == i) {
          // Distribution dist = model.getWeights(states[i]);
          Distribution dist = model.getWeights(states[i]);

          initG.add( ByteCode.make_aload           (modelLV));
          initG.add( ByteCode.make_aload           (statesLV));
          initG.add( ByteCode.make_iconst          (i));
          initG.add( ByteCode.make_aaload          ());
          initG.add( ByteCode.make_invokeinterface (_MarkovModel_getWeights));
          initG.add( ByteCode.make_astore          (distLV));

          int size = ((FiniteAlphabet) dist.getAlphabet()).size();
          Symbol[] transitionSymbols = new Symbol[size];

          initG.add( ByteCode.make_aload           (distLV));
          initG.add( ByteCode.make_invokeinterface (_Distribution_getAlphabet));
          initG.add( ByteCode.make_invokeinterface (_FiniteAlphabet_size));

          // t_i = new double[size]; // leave t_i on stack
          initG.add( ByteCode.make_iconst     (size));
          initG.add( ByteCode.make_newarray   (CodeUtils.TYPE_DOUBLE));
          initG.add( ByteCode.make_dup        ());
          initG.add( ByteCode.make_aload      (__init.getThis()));
          initG.add( ByteCode.make_swap       ());
          initG.add( ByteCode.make_putfield   (transitionFields[i]));

          // t_i[j] = scoreType.calculateScore(
          //            dist, states[ jj ]
          //          );
          // for each jj in dist - j is index from 0 in t_i in same order as jj
          int j = 0;
          for(int jj = 0; jj < states.length; jj++) {
            State state = states[jj];
            if(dist.getAlphabet().contains(state)) {
              transitionSymbols[j] = state;
              initG.add( ByteCode.make_dup             ());
              initG.add( ByteCode.make_iconst          (j));
              initG.add( ByteCode.make_aload           (__init.getVariable("scoreType")));
              initG.add( ByteCode.make_aload           (distLV));
              initG.add( ByteCode.make_aload           (statesLV));
              initG.add( ByteCode.make_iconst          (jj));
              initG.add( ByteCode.make_aaload          ());
              initG.add( ByteCode.make_invokeinterface (_ScoreType_calculateScore));
              initG.add( ByteCode.make_invokestatic    (_Math_log));
              initG.add( ByteCode.make_dastore         ());
              j++;
            }
          }

          transitionIndexers[i] = AlphabetManager.getAlphabetIndex(transitionSymbols);
          initG.add( ByteCode.make_pop ());
        }
      }

      // return nothing
      initG.add( ByteCode.make_return        ());

      return initG;
  }

  private InstructionVector createTransitionLastSum(
    GeneratedCodeMethod method,
    int i, int j,
    CodeField transition[], LocalVariable lastScore
  ) throws CodeException {
    InstructionVector sumV = new InstructionVector();

    // transition_j[i];
    sumV.add( ByteCode.make_aload    (method.getThis()));
    sumV.add( ByteCode.make_getfield (transition[j]));
    sumV.add( ByteCode.make_iconst   (i));
    sumV.add( ByteCode.make_daload   ());

    // lastScore[j]
    sumV.add( ByteCode.make_aload  (lastScore));
    sumV.add( ByteCode.make_iconst (j));
    sumV.add( ByteCode.make_daload  ());

    // add
    sumV.add( ByteCode.make_dadd   ());

    return sumV;
  }

  private CodeGenerator message(String message) {
    try {
      CodeClass _String = IntrospectedCodeClass.forClass(String.class);
      CodeClass _System = IntrospectedCodeClass.forClass(System.class);
      CodeField _System_out = _System.getFieldByName("out");
      CodeClass _PrintStream = IntrospectedCodeClass.forClass(PrintStream.class);
      CodeMethod _PrintStream_println = _PrintStream.getMethod("println", new CodeClass[] { _String });

      InstructionVector iv = new InstructionVector();

      iv.add( ByteCode.make_getstatic     (_System_out));
      iv.add( ByteCode.make_sconst        (message));
      iv.add( ByteCode.make_invokevirtual (_PrintStream_println));

      return iv;
    } catch (NoSuchFieldException nsfe) {
      throw new BioError( "Can't make message statements", nsfe);
    } catch (NoSuchMethodException nsme) {
      throw new BioError( "Can't make message statements", nsme);
    }
  }

  private CodeGenerator message(String message, LocalVariable var) {
    try {
      CodeClass _Object = IntrospectedCodeClass.forClass(Object.class);
      CodeClass _String = IntrospectedCodeClass.forClass(String.class);
      CodeClass _System = IntrospectedCodeClass.forClass(System.class);
      CodeField _System_out = _System.getFieldByName("out");
      CodeClass _PrintStream = IntrospectedCodeClass.forClass(PrintStream.class);
      CodeMethod _PrintStream_print = _PrintStream.getMethod("print", new CodeClass[] { _String });

      InstructionVector iv = new InstructionVector();

      iv.add( ByteCode.make_getstatic     (_System_out));
      iv.add( ByteCode.make_sconst        (message));
      iv.add( ByteCode.make_invokevirtual (_PrintStream_print));

      iv.add( ByteCode.make_getstatic     (_System_out));

      CodeMethod _PrintStream_printXln;
      if(var.getType().isPrimitive()) {
        _PrintStream_printXln = _PrintStream.getMethod("println", new CodeClass[] { var.getType() });
        if(var.getType() == CodeUtils.TYPE_INT) {
          iv.add( ByteCode.make_iload(var));
        } else if(var.getType() == CodeUtils.TYPE_DOUBLE) {
          iv.add( ByteCode.make_dload(var));
        } else {
          throw new BioError("Unsupported primative " + var.getType());
        }
      } else {
        iv.add( ByteCode.make_aload(var));
        _PrintStream_printXln = _PrintStream.getMethod("println", new CodeClass [] { _Object });
      }
      iv.add( ByteCode.make_invokevirtual (_PrintStream_printXln));

      return iv;
    } catch (NoSuchFieldException nsfe) {
      throw new BioError( "Can't make message statements", nsfe);
    } catch (NoSuchMethodException nsme) {
      throw new BioError( "Can't make message statements", nsme);
    } catch (CodeException ce) {
      throw new BioError( "Can't make message statements", ce);
    }
  }

  private CodeGenerator debug(CodeGenerator child) {
    if(debug) {
      return child;
    } else {
      return CodeUtils.DO_NOTHING;
    }
  }

  private int[] getAdvance(State s) {
    if(s instanceof EmissionState) {
      return ((EmissionState) s).getAdvance();
    } else {
      return new int[] {0, 0};
    }
  }

  private LocalVariable getAdvanced(LocalVariable[][] what, int[] advance) {
    return what[advance[0]][advance[1]];
  }

  private static class Factory implements CellCalculatorFactory {
    private final DP dp;
    private final Constructor forwards;
    private final Constructor backwards;
    private final Constructor viterbi;

    public Factory(
      DP dp,
      Constructor forwards,
      Constructor backwards,
      Constructor viterbi
    ) {
      this.dp = dp;
      this.viterbi = viterbi;
      this.forwards = forwards;
      this.backwards = backwards;
    }

    public CellCalculator forwards(ScoreType scoreType)
    throws
      IllegalSymbolException,
      IllegalAlphabetException,
      IllegalTransitionException
    {
      try {
        return (CellCalculator) forwards.newInstance(new Object[] { dp, scoreType });
      } catch (InstantiationException ie) {
        throw new BioError("Counld not instantiate auto-generated class");
      } catch (IllegalAccessException ie) {
        throw new BioError("Counld not instantiate auto-generated class");
      } catch (InvocationTargetException ie) {
        throw new BioError("Counld not instantiate auto-generated class");
      }
    }

    public CellCalculator backwards(ScoreType scoreType)
    throws
      IllegalSymbolException,
      IllegalAlphabetException,
      IllegalTransitionException
    {
      try {
        return (CellCalculator) backwards.newInstance(new Object[] { dp, scoreType });
      } catch (InstantiationException ie) {
        throw new BioError("Counld not instantiate auto-generated class");
      } catch (IllegalAccessException ie) {
        throw new BioError("Counld not instantiate auto-generated class");
      } catch (InvocationTargetException ie) {
        throw new BioError("Counld not instantiate auto-generated class");
      }
    }

    public CellCalculator viterbi(ScoreType scoreType, BackPointer terminal)
    throws
      IllegalSymbolException,
      IllegalAlphabetException,
      IllegalTransitionException
    {
      try {
        return (CellCalculator) viterbi.newInstance(new Object[] { dp, scoreType, terminal });
      } catch (InstantiationException ie) {
        throw new BioError( "Counld not instantiate auto-generated class", ie);
      } catch (IllegalAccessException ie) {
        throw new BioError( "Counld not instantiate auto-generated class-", ie);
      } catch (InvocationTargetException ie) {
        throw new BioError("Counld not instantiate auto-generated class using " + viterbi + " with " + dp + ", " + scoreType + ", " + terminal, ie);
      }
    }
  }
}
