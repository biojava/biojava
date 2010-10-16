

package org.biojava.utils.automata;

import java.util.Iterator;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;

public class PatternMaker
{
    // the NFA that will be the result of parsing the specified pattern.
    private Nfa nfa;
    private Tokenizer toke;
    private SymbolTokenization symtoke;

    private static class Tokenizer
    {
        private String packedTxt;
        private int ptr = 0;

        final static int EOL = -1;
        final static int SYMBOL_TOKEN = 0;
        final static int NUMERIC = 1;
        final static int LEFT_BRACE = 2;
        final static int RIGHT_BRACE = 3;
        final static int COMMA = 4;
        final static int LEFT_BRACKET = 5;
        final static int RIGHT_BRACKET = 6;
        final static int PLUS = 7;
        final static int ASTERISK = 8;
        final static int VERT_BAR = 9;
        final static int UNKNOWN = 999;

        private Tokenizer(String target)
        {
            packedTxt = pack(target);
        }

        private char getToken()
            throws IndexOutOfBoundsException
        {
            if (hasNext())
                 return packedTxt.charAt(ptr++);
            else
                throw new IndexOutOfBoundsException("text length: " + packedTxt.length() + " index: " + ptr);
        }

        private char peekToken()
            throws IndexOutOfBoundsException
        {
            if (hasNext())
                 return packedTxt.charAt(ptr);
            else
                throw new IndexOutOfBoundsException("text length: " + packedTxt.length() + " index: " + ptr);
        }

        private int nextTokenType()
        {
            if (!hasNext()) return EOL;

            char nextCh = packedTxt.charAt(ptr);

            // symbol tokens assumed to be alphas.
            if (Character.isLetter(nextCh))
                return SYMBOL_TOKEN;
            if (Character.isDigit(nextCh))
                return NUMERIC;
            // now check for specific chars
            if (nextCh == '{')
                return LEFT_BRACE;
            if (nextCh == '}')
                return RIGHT_BRACE;
            if (nextCh == ',')
                return COMMA;
            if (nextCh == '(')
                return LEFT_BRACKET;
            if (nextCh == ')')
                return RIGHT_BRACKET;
            if (nextCh == '+')
                return PLUS;
            if (nextCh == '*')
                return ASTERISK;
            if (nextCh == '|')
                return VERT_BAR;
            return UNKNOWN;
        }

        private boolean hasNext()
        {
            return ptr < packedTxt.length();
        }

        /**
         * produces a version of the String with whitespace removed.
         */
        private String pack(String source)
        {
            StringBuffer packedString = new StringBuffer();

            for (int i=0; i < source.length(); i++) {
                char currCh;
                if (!Character.isWhitespace(currCh = source.charAt(i))) {
                    packedString.append(currCh);
                }
            }

            return packedString.toString();
        }
    }

    private static class Range
    {
        private int min = 1;
        private int max = 1;
        private int getMin() { return min; }
        private int getMax() { return max; }

        private Range(int min, int max)
        {
            this.min = min;
            this.max = max;
        }

        private boolean once() { return (min == 1) && (max == 1); }
    }

    PatternMaker(String patternString, FiniteAlphabet alfa)
        throws BioException
    {
        toke = new Tokenizer(patternString);
        nfa = new Nfa(patternString, alfa);
        symtoke = alfa.getTokenization("token");
    }

    /**
     * Compiles the regex described by patternString into a DFA.
     * The DFA can subsequently be converted into a state
     * state machine for searching. 
     * <b>This is the main work method for this class.</b>
     */
    static FiniteAutomaton compilePattern(String patternString, FiniteAlphabet alfa)
        throws BioException, AutomatonException
    {
        PatternMaker maker = new PatternMaker(patternString, alfa);
        return maker.parse();
    }

    private FiniteAutomaton parse()
        throws AutomatonException
    {
        NfaSubModel result = parse(nfa);

        // complete linking up model to the start/end.
        nfa.addEpsilonTransition(nfa.getStart(), result.getStart());
        nfa.addEpsilonTransition(result.getEnd(), nfa.getEnd());
        //System.out.println(nfa.toString());
        nfa.doEpsilonClosure();
        //System.out.println(nfa.toString());
        return new DfaBuilder(nfa).getDFA();
    }    

    private NfaSubModel parse(NfaBuilder delegate)
        throws AutomatonException
    {
        /**
         * The model is that the pattern here can have
         * alternative patterns so a single pattern ends
         * up being a pattern of branching 1.
         */
        NfaSubModel returnSubModel = new NfaSubModel(delegate);
        NfaSubModel branchSubModel = new NfaSubModel(returnSubModel);
        // make it link silently so append works properly.
        branchSubModel.addEpsilonTransition(branchSubModel.getStart(), branchSubModel.getEnd());
        // do I have content in this branch?
        boolean gotContent =false;

        Range times;
        while (toke.hasNext()) {
            int tokenType = toke.nextTokenType();

            NfaSubModel currSubModel = null;

            switch (tokenType) {
                case Tokenizer.LEFT_BRACKET:
                    //System.out.println("processing left bracket" + toke.peekToken());
                    gotContent = true;
                    currSubModel = parseSubPattern(branchSubModel);
                    times = parseIterations();
                    currSubModel = reiterate(currSubModel, branchSubModel,  times);
                    branchSubModel.append(currSubModel);
                    break;
                case Tokenizer.SYMBOL_TOKEN:
                    //System.out.println("processing symbol " + toke.peekToken());
                    gotContent = true;
                    currSubModel = parseSymbol(branchSubModel);
                    times = parseIterations();
                    currSubModel = reiterate(currSubModel, branchSubModel, times);
//                    System.out.println("before\n" + branchSubModel);
                    branchSubModel.append(currSubModel);
//                    System.out.println("after\n" + branchSubModel);
                    break;
                case Tokenizer.VERT_BAR:
                    //System.out.println("processing bar" + toke.peekToken());
                    // link current branch into return value
                    if (!gotContent) throw new AutomatonException("no content in this branch!");
                    System.out.println(returnSubModel.getStart().getID() + " " + branchSubModel.getStart().getID());
                    returnSubModel.addEpsilonTransition(
                        returnSubModel.getStart(),
                        branchSubModel.getStart());
                    returnSubModel.addEpsilonTransition(
                        branchSubModel.getEnd(),
                        returnSubModel.getEnd());
                    // start new branch
                    branchSubModel = new NfaSubModel(returnSubModel);
                    branchSubModel.addEpsilonTransition(branchSubModel.getStart(), branchSubModel.getEnd());
                    gotContent = false;
                    toke.getToken();
                    break;
                case Tokenizer.RIGHT_BRACKET:
                    //System.out.println("processing right bracket" + toke.peekToken());
                    // link current branch into return value
                    if (!gotContent) throw new AutomatonException("no content in this branch!");
                    //System.out.println(returnSubModel.getStart().getID() + " " + branchSubModel.getStart().getID());
                    returnSubModel.addEpsilonTransition(
                        returnSubModel.getStart(),
                        branchSubModel.getStart());
                    returnSubModel.addEpsilonTransition(
                        branchSubModel.getEnd(),
                        returnSubModel.getEnd());
                    return returnSubModel; // note that the right bracket is consumed by caller.
                default:
                    throw new AutomatonException("Illegal symbol encountered: " + toke.peekToken());
            }
        }

        if (gotContent) {
            // the latest content would be in branchSubModel
            // and this needs to be connected back into the
            // returnSubModel.
            returnSubModel.addEpsilonTransition(
                returnSubModel.getStart(),
                branchSubModel.getStart());
            returnSubModel.addEpsilonTransition(
                branchSubModel.getEnd(),
                returnSubModel.getEnd());
            return returnSubModel;
        }
        else throw new AutomatonException("no content in this branch!");
    }

    NfaSubModel parseSubPattern(NfaBuilder delegate)
        throws AutomatonException
    {
        NfaSubModel returnSubModel;

        // consume the left bracket
        toke.getToken();
        returnSubModel = parse(delegate);
        // consume the right bracket
        if (toke.nextTokenType() != Tokenizer.RIGHT_BRACKET)
            throw new AutomatonException("Missing right bracket!");
        else
            toke.getToken();

        return returnSubModel;
    }

    NfaSubModel parseSymbol(NfaBuilder delegate)
        throws AutomatonException
    {
        NfaSubModel returnSubModel = new NfaSubModel(delegate);

        Symbol sym;
        try {
            sym = symtoke.parseToken(Character.toString(toke.getToken()));
        }
        catch (IllegalSymbolException ise) {
            throw new AutomatonException(ise);
        }

        FiniteAutomaton.Node pre = returnSubModel.addNode(false);
        FiniteAutomaton.Node post = returnSubModel.addNode(false);
        if (sym instanceof AtomicSymbol) {
            returnSubModel.addTransition(pre, post, sym);
            returnSubModel.addEpsilonTransition(returnSubModel.getStart(), pre);
            returnSubModel.addEpsilonTransition(post, returnSubModel.getEnd());
        }
        else {
            for (Iterator symI = ((FiniteAlphabet) sym.getMatches()).iterator();
                symI.hasNext(); ) {
                Symbol atomicSym = (Symbol) symI.next();

                returnSubModel.addTransition(pre, post, atomicSym);
                returnSubModel.addEpsilonTransition(returnSubModel.getStart(), pre);
                returnSubModel.addEpsilonTransition(post, returnSubModel.getEnd());
            }
        }

        return returnSubModel;
    }

    Range parseIterations()
        throws AutomatonException
    {
        int tokenType = toke.nextTokenType();

        switch (tokenType) {
            case Tokenizer.LEFT_BRACE:
                return getIterations();
            case Tokenizer.PLUS:
                toke.getToken();
                return new Range(1, Integer.MAX_VALUE);
            case Tokenizer.ASTERISK:
                toke.getToken();
                return new Range(0, Integer.MAX_VALUE);
            default:
                return new Range(1, 1);
        }

    }

    private Range getIterations()
        throws AutomatonException
    {
        int min = 0;
        int max = 0;

        // consume the left brace
        toke.getToken();

        // there can either be one or two numbers
        boolean onSecondArg = false;
        StringBuffer numString = new StringBuffer();

        while (toke.hasNext()) {
            int tokenType = toke.nextTokenType();

            switch (tokenType) {
                case Tokenizer.NUMERIC:
                    //System.out.println("adding symbol");
                    if (numString == null) numString = new StringBuffer();
                    numString.append(toke.getToken());
                    break;
                case Tokenizer.COMMA:
                    toke.getToken();
                    if (!onSecondArg) {
                        //System.out.println("numString is " + numString);
                        if (numString.length() > 0)
                            min = Integer.parseInt(numString.toString());
                        else
                            min = 0;
                        numString = null;
                        onSecondArg = true;
                    }
                    else {
                        throw new AutomatonException("only two arguments permitted.");
                    }
                    break;
                case Tokenizer.RIGHT_BRACE:
                    toke.getToken();
                    if (onSecondArg) {
                        if (numString != null)
                            max = Integer.parseInt(numString.toString());
                        else
                            max = Integer.MAX_VALUE;
                    }
                    else {
                        min = max = Integer.parseInt(numString.toString());
                    }
                    return new Range(min, max);
                default:
                    throw new AutomatonException(toke.getToken() + " is not valid in an iteration specifier.");
            }
        }

        throw new AutomatonException("unexpected error.");
    }


    NfaSubModel reiterate(NfaSubModel base, NfaBuilder delegate, Range count)
    {
        // no repeats necessary!
        if (count.once()) return base;

        NfaSubModel returnSubModel = new NfaSubModel(delegate);
        FiniteAutomaton.Node lastNode = returnSubModel.getStart();

        // fulfill min number of reiterations
        for (int i=0; i < count.getMin(); i++) {
            NfaSubModel duplSubModel = base.duplicate();
            returnSubModel.addEpsilonTransition(lastNode, duplSubModel.getStart());
            lastNode = duplSubModel.getEnd();
        }

        // is there a finite upper bound?
        if (count.getMax() != Integer.MAX_VALUE) {
            for (int i = count.getMin(); i < count.getMax(); i++) {
                returnSubModel.addLambdaTransition(lastNode, returnSubModel.getEnd());
                NfaSubModel duplSubModel = base.duplicate();
                returnSubModel.addEpsilonTransition(lastNode, duplSubModel.getStart());
                lastNode = duplSubModel.getEnd();
            }
        }
        else {
            // infinite upper bound
            NfaSubModel duplSubModel = base.duplicate();
            returnSubModel.addEpsilonTransition(lastNode, duplSubModel.getStart());
            returnSubModel.addEpsilonTransition(duplSubModel.getEnd(), lastNode);
        }

        returnSubModel.addEpsilonTransition(lastNode, returnSubModel.getEnd());

        return returnSubModel;
    }
}
