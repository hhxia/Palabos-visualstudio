// LBMDoc.h : interface of the CLBMDoc class
//


#pragma once


class CLBMDoc : public CDocument
{
protected: // create from serialization only
	CLBMDoc();
	DECLARE_DYNCREATE(CLBMDoc)

// Attributes
public:

// Operations
public:

// Overrides
public:
	virtual BOOL OnNewDocument();
	virtual void Serialize(CArchive& ar);

// Implementation
public:
	virtual ~CLBMDoc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// Generated message map functions
protected:
	DECLARE_MESSAGE_MAP()
};


