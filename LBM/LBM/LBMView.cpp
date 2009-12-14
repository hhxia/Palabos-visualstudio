// LBMView.cpp : implementation of the CLBMView class
//

#include "stdafx.h"
#include "LBM.h"

#include "LBMDoc.h"
#include "LBMView.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CLBMView

IMPLEMENT_DYNCREATE(CLBMView, CView)

BEGIN_MESSAGE_MAP(CLBMView, CView)
	// Standard printing commands
	ON_COMMAND(ID_FILE_PRINT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_DIRECT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_PREVIEW, &CView::OnFilePrintPreview)
END_MESSAGE_MAP()

// CLBMView construction/destruction

CLBMView::CLBMView()
{
	// TODO: add construction code here

}

CLBMView::~CLBMView()
{
}

BOOL CLBMView::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: Modify the Window class or styles here by modifying
	//  the CREATESTRUCT cs

	return CView::PreCreateWindow(cs);
}

// CLBMView drawing

void CLBMView::OnDraw(CDC* /*pDC*/)
{
	CLBMDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	if (!pDoc)
		return;

	// TODO: add draw code for native data here
}


// CLBMView printing

BOOL CLBMView::OnPreparePrinting(CPrintInfo* pInfo)
{
	// default preparation
	return DoPreparePrinting(pInfo);
}

void CLBMView::OnBeginPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: add extra initialization before printing
}

void CLBMView::OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: add cleanup after printing
}


// CLBMView diagnostics

#ifdef _DEBUG
void CLBMView::AssertValid() const
{
	CView::AssertValid();
}

void CLBMView::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CLBMDoc* CLBMView::GetDocument() const // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CLBMDoc)));
	return (CLBMDoc*)m_pDocument;
}
#endif //_DEBUG


// CLBMView message handlers
