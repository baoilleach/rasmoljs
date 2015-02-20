/* tokens.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, December 1998
 * Version 2.6.4
 */

/* Lexeme Tokens */
#define IdentTok       256
#define NumberTok      257
#define FloatTok       258
#define StringTok      259

/* Command Tokens */
#define AdviseTok      260
#define BackboneTok    261
#define CartoonTok     262
#define CentreTok      263
#define ClipboardTok   264
#define ColourTok      265
#define ConnectTok     266
#define DashTok        267
#define DefineTok      268
#define DelayTok       269
#define DisplayTok     270
#define EchoTok        271
#define ExitTok        272
#define HelpTok        273
#define LabelTok       274
#define LoadTok        275
#define LoopTok        276
#define MonitorTok     277
#define MoveTok        278
#define PrintTok       279
#define QuitTok        280
#define RefreshTok     281
#define RenumTok       282
#define ResetTok       283
#define ResizeTok      284
#define RestoreTok     285
#define RestrictTok    286
#define RotateTok      287
#define SaveTok        288
#define ScriptTok      289
#define SelectTok      290
#define SetTok         291
#define ShowTok        292
#define SlabTok        293
#define SourceTok      294
#define SpacefillTok   295
#define StructureTok   296
#define SymmetryTok    297
#define TitleTok       298
#define TraceTok       299
#define TranslateTok   300
#define ViewTok        301
#define WaitTok        302
#define WireframeTok   303
#define WriteTok       304
#define ZapTok         305
#define ZoomTok        306

/* Predicate Tokens */
#define IsPredTok(x)   (((x)>=310) && ((x)<=348))
#define PredTokOrd(x)  ((x)-310)
#define PredTokChr(x)  ((x)+310)

#define AlphaTok       310
#define AminoTok       311
#define ATTok          312
#define BondedTok      313
#define CGTok          314
#define CystineTok     315
#define DNATok         316
#define HelixTok       317
#define HeteroTok      318
#define HydrogenTok    319
#define IonTok         320
#define LigandTok      321
#define MainChainTok   322
#define NucleicTok     323
#define ProteinTok     324
#define PurineTok      325
#define PyrimidineTok  326
#define RNATok         327
#define SelectedTok    328
#define SheetTok       329
#define SidechainTok   330
#define SolventTok     331
#define TurnTok        332
#define WaterTok       333

#define AcidicTok      334
#define AcyclicTok     335
#define AliphaticTok   336
#define AromaticTok    337
#define BasicTok       338
#define BuriedTok      339
#define ChargedTok     340
#define CyclicTok      341
#define HydrophobicTok 342
#define LargeTok       343
#define MediumTok      344
#define NeutralTok     345
#define PolarTok       346
#define SmallTok       347
#define SurfaceTok     348


/* Property Tokens */
#define IsPropTok(x)   (((x)>=350) && ((x)<=355))
#define TemperatureTok 350
#define RadiusTok      351
#define AtomNoTok      352
#define ElemNoTok      353
#define ModelTok       354
#define ResNoTok       355

/* File Format Tokens */
/* Warning! Tokens are related to Format values */
#define IsMoleculeToken(x)  (((x)>=360) && ((x)<=375))

#define PDBTok         360
#define MacroModelTok  361
#define GaussianTok    362
#define AlchemyTok     363
#define NMRPDBTok      364
#define CharmmTok      365
#define BiosymTok      366
#define MOPACTok       367
#define SHELXTok       368
#define Mol2Tok        369
#define FDATTok        370
#define MMDBTok        371
#define MDLTok         372
#define XYZTok         373
#define CIFTok         374
#define CEXTok         375

#define GraspTok       376

/* Raster Tokens */
#define IsImageToken(x) (((x)>=380) && ((x)<=398))
#define GIFTok         380
#define PPMTok         381
#define SUNTok         382
#define SUNRLETok      383
#define EPSFTok        384
#define PICTTok        385
#define IRISTok        386
#define BMPTok         387
#define MonoPSTok      388
#define JPEGTok        389
#define PNGTok         390

#define VectPSTok      391
#define KinemageTok    392
#define MolScriptTok   393
#define POVRayTok      394
#define VRMLTok        395
#define STLTok         396
#define DXFTok         397
#define PLYTok         398

/* Feature Tokens */
#define AtomTok        400
#define BondTok        401
#define DotsTok        402
#define HBondTok       403
#define RibbonTok      404
#define SSBondTok      405
#define Ribbon1Tok     406
#define Ribbon2Tok     407

/* Expression Tokens */
#define TrueTok        410
#define FalseTok       411
#define AllTok         412
#define NoneTok        413
#define AndTok         414
#define OrTok          415
#define NotTok         416
#define WithinTok      417
#define XorTok         418

/* Colour Tokens */
/* Warning! Tokens are related to colour values */
#define IsColourToken(x) (((x)>=420) && ((x)<=443))
#define Token2Colour(x)  ((x)-420)

#define BlackTok       420
#define BlueTok        421
#define BlueTintTok    422
#define BrownTok       423
#define CyanTok        424
#define GoldTok        425
#define GrayTok        426
#define GreenTok       427
#define GreenBlueTok   428
#define GreenTintTok   429
#define HotPinkTok     430
#define MagentaTok     431
#define OrangeTok      432
#define PinkTok        433
#define PinkTintTok    434
#define PurpleTok      435
#define RedTok         436
#define RedOrangeTok   437
#define SeaGreenTok    438
#define SkyBlueTok     439
#define VioletTok      440
#define WhiteTok       441
#define YellowTok      442
#define YellowTintTok  443

#define CPKTok         444
#define ShapelyTok     445
#define ResidueTok     446
#define UserTok        447
#define GroupTok       448
#define ChainTok       449
#define TypeTok        440
#define PotentialTok   451
#define ChargeTok      452

/* Variable Tokens */
#define AmbientTok     460
#define AxesTok        461
#define BackFadeTok    462
#define BackgroundTok  463
#define BondModeTok    464
#define BoundBoxTok    465
#define DepthCueTok    466
#define FontSizeTok    467
#define HourGlassTok   468
#define MenusTok       469
#define MouseTok       470
#define PickingTok     471
#define ShadowTok      472
#define SlabModeTok    473
#define SpecularTok    474
#define SpecPowerTok   475
#define StrandsTok     476
#define TransparentTok 477
#define UnitCellTok    478

/* SlabMode Tokens */
#define RejectTok      480
#define HalfTok        481
#define HollowTok      482
#define SolidTok       483
#define SectionTok     484

/* MouseMode Tokens */
#define RasMolTok      485
#define InsightTok     486
#define QuantaTok      487
#define SybylTok       488

/* Information Tokens */
#define InfoTok        490
#define SequenceTok    491
#define VersionTok     492

/* Display Mode Tokens */
#define NormalTok      495
#define StereoTok      496
#define MonoTok        497
#define HardwareTok    498

/* Axis Tokens */
#define XTok           500
#define YTok           501
#define ZTok           502

/* Picking Tokens */
#define IdentifyTok    505
#define DistanceTok    506
#define AngleTok       507
#define TorsionTok     508
#define OriginTok      509

/* Resolution Tokens */
#define HighTok        510
#define LowTok         511

/* Surface Tokens */
#define LineTok        512
#define VectorTok      513
#define TranslucentTok 514

/* Misc Tokens */
#define InLineTok      515
#define VDWTok         516

int LookUpKeyword( char *ptr );

