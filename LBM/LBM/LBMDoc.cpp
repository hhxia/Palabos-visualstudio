// LBMDoc.cpp : implementation of the CLBMDoc class
//

#include "stdafx.h"
#include "LBM.h"

#include "LBMDoc.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CLBMDoc

IMPLEMENT_DYNCREATE(CLBMDoc, CDocument)

BEGIN_MESSAGE_MAP(CLBMDoc, CDocument)
END_MESSAGE_MAP()


// CLBMDoc construction/destruction

CLBMDoc::CLBMDoc()
{
	// TODO: add one-time construction code here

}

CLBMDoc::~CLBMDoc()
{
}

BOOL CLBMDoc::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;

	// TODO: add reinitialization code here
	// (SDI documents will reuse this document)

	return TRUE;
}




// CLBMDoc serialization

void CLBMDoc::Serialize(CArchive& ar)
{
	if (ar.IsStoring())
	{
		// TODO: add storing code here
	}
	else
	{
		// TODO: add loading code here
	}
}


// CLBMDoc diagnostics

#ifdef _DEBUG
void CLBMDoc::AssertValid() const
{
	CDocument::AssertValid();
}

void CLBMDoc::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG


// CLBMDoc commands
