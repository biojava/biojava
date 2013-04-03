// $Id: Property.java,v 1.17 2009/10/26 23:29:39 cmzmasek Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// Copyright (C) 2000-2001 Washington University School of Medicine
// and Howard Hughes Medical Institute
// All rights reserved
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: cmzmasek@yahoo.com
// WWW: www.phylosoft.org/forester

package org.forester.phylogeny.data;

import java.io.IOException;
import java.io.Writer;
import java.util.StringTokenizer;

import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.util.ForesterUtil;

public class Property implements PhylogenyData {

    private final String    _value;
    private final String    _ref;
    private final String    _unit;
    private final String    _datatype;
    private final AppliesTo _applies_to;
    private final String    _id_ref;

    public Property( final String ref,
                     final String value,
                     final String unit,
                     final String datatype,
                     final AppliesTo applies_to ) {
        this( ref, value, unit, datatype, applies_to, "" );
    }

    private Property( final String ref,
                      final String value,
                      final String unit,
                      final String datatype,
                      final AppliesTo applies_to,
                      final boolean dummy ) {
        _ref = ref;
        _value = value;
        _unit = unit;
        _datatype = datatype;
        _applies_to = applies_to;
        _id_ref = "";
    }

    public Property( final String ref,
                     final String value,
                     final String unit,
                     final String datatype,
                     final AppliesTo applies_to,
                     final String id_ref ) {
        if ( !ForesterUtil.isEmpty( ref ) && ( ref.indexOf( ":" ) < 1 ) ) {
            throw new IllegalArgumentException( "property reference [" + ref
                    + "] is not in the expected format (missing a \":\")" );
        }
        if ( !ForesterUtil.isEmpty( unit ) && ( unit.indexOf( ":" ) < 1 ) ) {
            throw new IllegalArgumentException( "property unit [" + unit
                    + "] is not in the expected format (missing a \":\")" );
        }
        if ( !ForesterUtil.isEmpty( datatype ) && ( datatype.indexOf( ":" ) < 1 ) ) {
            throw new IllegalArgumentException( "property datatype [" + unit
                    + "] is not in the expected format (missing a \":\")" );
        }
        _ref = ref;
        _value = value;
        _unit = unit;
        _datatype = datatype;
        _applies_to = applies_to;
        _id_ref = id_ref;
    }

    @Override
    public StringBuffer asSimpleText() {
        return new StringBuffer( getValue() );
    }

    @Override
    public StringBuffer asText() {
        final StringBuffer sb = new StringBuffer();
        sb.append( getRef() );
        sb.append( ": " );
        sb.append( getValue() );
        if ( !ForesterUtil.isEmpty( getUnit() ) ) {
            sb.append( getUnit() );
        }
        return sb;
    }

    @Override
    public PhylogenyData copy() {
        return new Property( new String( getRef() ),
                             new String( getValue() ),
                             new String( getUnit() ),
                             new String( getDataType() ),
                             getAppliesTo(),
                             new String( getIdRef() ) );
    }

    public AppliesTo getAppliesTo() {
        return _applies_to;
    }

    public String getDataType() {
        return _datatype;
    }

    public String getIdRef() {
        return _id_ref;
    }

    public String getRef() {
        return _ref;
    }

    public String getUnit() {
        return _unit;
    }

    public String getValue() {
        return _value;
    }

    @Override
    public boolean isEqual( final PhylogenyData data ) {
        if ( data == null ) {
            return false;
        }
        return ( ( Property ) data ).getValue().equals( getValue() )
                && ( ( Property ) data ).getUnit().equals( getUnit() )
                && ( ( Property ) data ).getRef().equals( getRef() );
    }

    @Override
    public StringBuffer toNHX() {
        final StringBuffer nhx = new StringBuffer();
        nhx.append( ":X" );
        switch ( getAppliesTo() ) {
            case CLADE:
                nhx.append( "C=" );
                break;
            case NODE:
                nhx.append( "N=" );
                break;
            case PARENT_BRANCH:
                nhx.append( "B=" );
                break;
            case PHYLOGENY:
                nhx.append( "P=" );
                break;
            case ANNOTATION:
                nhx.append( "S=" );
                break;
            default:
                nhx.append( "O=" );
                break;
        }
        if ( !getDataType().equals( "" ) ) {
            if ( getDataType().equals( "xsd:string" ) ) {
                nhx.append( "S=" );
            }
            else if ( getDataType().equals( "xsd:long" ) ) {
                nhx.append( "L=" );
            }
            else if ( getDataType().equals( "xsd:decimal" ) ) {
                nhx.append( "D=" );
            }
            else if ( getDataType().equals( "xsd:boolean" ) ) {
                nhx.append( "B=" );
            }
            else if ( getDataType().equals( "xsd:anyUR" ) ) {
                nhx.append( "U=" );
            }
        }
        nhx.append( getRef() );
        nhx.append( "=" );
        nhx.append( getValue() );
        if ( !getUnit().equals( "" ) ) {
            nhx.append( "=" );
            nhx.append( getUnit() );
        }
        return nhx;
    }

    @Override
    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        PhylogenyDataUtil.appendElement( writer,
                                         PhyloXmlMapping.PROPERTY,
                                         getValue(),
                                         PhyloXmlMapping.PROPERTY_REF,
                                         getRef(),
                                         PhyloXmlMapping.PROPERTY_UNIT,
                                         getUnit(),
                                         PhyloXmlMapping.PROPERTY_DATATYPE,
                                         getDataType(),
                                         PhyloXmlMapping.PROPERTY_APPLIES_TO,
                                         getAppliesTo().toString(),
                                         PhyloXmlMapping.ID_REF,
                                         getIdRef(),
                                         indentation );
    }

    @Override
    public String toString() {
        return asText().toString();
    }

    public static Property createFromNhxString( final String nhx ) throws IllegalArgumentException {
        final StringTokenizer st = new StringTokenizer( nhx, "=" );
        final int tokens = st.countTokens();
        final String error = "error in NHX property tag format: "
                + "expected: X[N|B|C|S|T|P|O]=<datatype>=<ref>=<value>[=<unit>], got: \"" + nhx + "\" instead";
        if ( ( tokens != 4 ) && ( tokens != 5 ) ) {
            throw new IllegalArgumentException( error );
        }
        final String first = st.nextToken();
        AppliesTo applies_to = null;
        if ( first.equals( "XN" ) ) {
            applies_to = AppliesTo.NODE;
        }
        else if ( first.equals( "XB" ) ) {
            applies_to = AppliesTo.PARENT_BRANCH;
        }
        else if ( first.equals( "XC" ) ) {
            applies_to = AppliesTo.CLADE;
        }
        else if ( first.equals( "XS" ) ) {
            applies_to = AppliesTo.ANNOTATION;
        }
        else if ( first.equals( "XT" ) ) {
            applies_to = AppliesTo.OTHER;
        }
        else if ( first.equals( "XP" ) ) {
            applies_to = AppliesTo.PHYLOGENY;
        }
        else if ( first.equals( "XO" ) ) {
            applies_to = AppliesTo.OTHER;
        }
        else {
            throw new IllegalArgumentException( error );
        }
        String datatype = st.nextToken();
        if ( datatype.equals( "S" ) ) {
            datatype = "xsd:string";
        }
        else if ( datatype.equals( "L" ) ) {
            datatype = "xsd:long";
        }
        else if ( datatype.equals( "D" ) ) {
            datatype = "xsd:decimal";
        }
        else if ( datatype.equals( "B" ) ) {
            datatype = "xsd:boolean";
        }
        else if ( datatype.equals( "U" ) ) {
            datatype = "xsd:anyURI";
        }
        final String ref = st.nextToken();
        final String value = st.nextToken();
        String unit = "";
        if ( tokens == 5 ) {
            unit = st.nextToken();
        }
        return new Property( ref, value, unit, datatype, applies_to, true );
    }

    public static enum AppliesTo {
        PHYLOGENY {

            @Override
            public String toString() {
                return "phylogeny";
            }
        },
        CLADE {

            @Override
            public String toString() {
                return "clade";
            }
        },
        NODE {

            @Override
            public String toString() {
                return "node";
            }
        },
        ANNOTATION {

            @Override
            public String toString() {
                return "annotation";
            }
        },
        PARENT_BRANCH {

            @Override
            public String toString() {
                return "parent_branch";
            }
        },
        OTHER {

            @Override
            public String toString() {
                return "other";
            }
        }
    }
}
